clear
sourcePath = 'D:\MGH DATA\SharingData';
addpath('D:\Single input intravascular polarimetry with mirror state\single input processing toolbox');
%%
% These are the pullbacks from which Ken selected the slices for his
% analysis.
Files = importdata('FilesNameForComparison.mat');
logLim = [55,110];
N = 5;
st = struct;
makeOth = 0;  % make orthogonal or not, Not by default
%% folder name define
if makeOth
    folderName = 'MP orthogonal';
else
    folderName = 'MP normal';
end
presentDir = pwd;
cd('D:\SingleInputProcessingResults\Mirror state')
mkdir(folderName)
cd(folderName)
st.fwx = 10;
st.makeOth = makeOth;
set(0,'DefaultFigureVisible', 'off');
for ind = 1:numel(Files)
    fInd = ind
    logF = readLogFile(fullfile(sourcePath,Files{fInd}(1:6),Files{fInd}));
    path = fullfile(sourcePath,Files{fInd}(1:6),Files{fInd});
    %% first let's generate the mirror state of the full spectrum
    S1S = [];S2S = [];
    for frameIdx = 1:5:50
        [S1,S2] = recstrTom(path,[1,1]*(frameIdx*5),st);
        int = log10(tom2Int(S1,S2));
        S1S = cat(2,permute(S1,[1 2 4 3]),S1S); S2S = cat(2,permute(S2,[1 2 4 3]),S2S);
    end
    [mp1f, mp2f, conts] = mirrorPointExtract(S1S, S2S, st);
    saveas(gcf,['full',num2str(fInd),'.bmp']);close; 
    MP1f(fInd,:) = mp1f; MP2f(fInd,:) = mp2f; Conf(fInd,:) = conts;
    %% Second, let's genearte the mirror state for spectral binning
    pst = st; pst.MP1 = mp1f; pst.MP2 = mp2f; % use the mirror state from full spectra as the see mirror state
    S1BS = []; S2BS = [];
    for frameIdx = 1:5:50
        % spectral binning stokes reconstruction
        st2 = st;
        st2.window = N;
        st2.skipLastPoints = 30;
        [S1,S2] = recstrTom(path,[1,1]*(frameIdx*5),st2);
        S1BS = cat(2,S1,S1BS); S2BS = cat(2,S2,S2BS);
    end
    [MP1, MP2, cond] = mirrorPointExtract(S1BS, S2BS, pst);
    Cons(fInd,:,:) = mat2gray(cond,[50 250]);
    saveas(gcf,['Binned',num2str(fInd),'.bmp']);close;
    MP1s(fInd,:,:) = MP1; MP2s(fInd,:,:) = MP2;
end
cd(presentDir);
save('FullSpectral Mirror state', 'MP1f', 'MP2f', 'Conf')
save('Binned Mirror state', 'MP1s', 'MP2s', 'Cons')
