clear
cd('D:\mirror state processing')
load('FilesR1.mat');
sourcePath = 'D:\MGH DATA\SharingData';
DataSetName = importdata('FilesNameForComparison.mat');
addpath('D:\mirror state processing\single input processing toolbox')
%
%
%% ===============================setting parameters for the processing===================================================================== 
%--choose what reliablility metric to be used
relMetric = 'sina';  % sina means using reliability metric associated with sina, and dist means using metric associated with distance metric
if strcmp(relMetric,'sina')
    load('D:\mirror state processing\Reliability metric generation\RelMetsina.mat')
    pst.relm = 'sina';
else
    load('D:\mirror state processing\Reliability metric generation\RelMetdistance.mat')
    pst.relm = 'dist';
end
pst.Apval = Apval; pst.Bpval = Bpval;
%chose approximate or accurate mode to reconstruct local reatardation
pst.mode = 'approximate'; %% approximate or accurate mode
%--make orthogonal or not, this will be removed afterwards
makeOrth = 0; %% no making orthognal.
%% =========================================================================================================================================
%
%
%%load mirror state
load('Binned Mirror state.mat')
outcomeName = pst.mode;
%%create a folder for the results repository
resultFolderName = ['Intravascular ', datestr(now,29), relMetric];
preDir = pwd; cd('D:\SingleInputProcessingResults');mkdir(resultFolderName);cd(resultFolderName);mkdir('Numerical');mkdir('images');cd(preDir);
%% -------------------------------------------------------------------------------------------------------------------
%%processing parameters setting
logLim = [55,110];
N = 5;
st = struct;
pst.fwx = 6;
pst.dz = 5;
pst.makeOrth = makeOrth;
%%Gaussian fitering kernel
h = filterGen(pst.fwx).';
Ret1 = []; Ret2 = []; Ret3 = [];Rel2 = []; Rel3 = [];
set(0,'DefaultFigureVisible', 'off');
MPindex = 1:43;
fileNo = 0;
for fInd = 1:5
    fInd
    path = fullfile(sourcePath,MsmtR{fInd}(1:6),MsmtR{fInd});
    temp = contains(DataSetName,string(MsmtR{fInd}));
    MPfIndex = MPindex(temp);
    %% mirror state
    pst.MP1 = squeeze(MP1s(MPfIndex,:,:));pst.MP2 = squeeze(MP2s(MPfIndex,:,:));
    [S1,S2] = recstrTom(path,[1,1]*SliceR(fInd),st);
    %%----filtering
    S1f = imfilter(S1,h,'circular');S2f = imfilter(S2,h,'circular');
    %%--normalization
    If1 = sqrt(dot(S1f,S1f,3));If2 = sqrt(dot(S2f,S2f,3));
    S1n = S1f./repmat(If1,[1,1,3]);S2n = S2f./repmat(If2,[1,1,3]);
    %%----detect the catheter------------------------
    int = tom2Int(S1,S2); cath = findCatheter(int);
    clipLimit = max(cath(3,:)); pst.cath = cath;
    % spectral binning reconstruction
    st2 = st; st2.window = N;st2.skipLastPoints = 30;pst.clipLimit = clipLimit;pst.offset = offsetR(fInd);pst.fwx = 6;
    [S1,S2] = recstrTom(path,[1,1]*SliceR(fInd),st2);
    %----- benchmark two input processing
    out2input = PSProcessLocal(S1,S2,pst);
    %------single input processing
    out1input = MirrorStateProcessAcc(S1,S2,pst);
%     out1input = MirrorStateProcess_NewSina(S1,S2,pst);
    %-------cmparison validataion
    out2 = comVisualization(int,out2input,out1input,pst,1);
    dirOig = pwd;
    cd(['D:\SingleInputProcessingResults\', resultFolderName]);
    cd('Numerical'); 
%     mkdir(num2str(fInd)); cd(num2str(fInd)); saveimages;cd('..\');
    saveas(gcf,[outcomeName,num2str(fInd),'.bmp']);
    close all; cd('..\'); cd('images')
    imwrite(cat(2,repmat(out2.cint,[1 1 3]),out2.cim1m, out2.cim2m, out2.cim3m),[outcomeName,num2str(fInd),'.bmp'],'BMP')
    cd(dirOig);
    Slope(fInd,:) = out2.slope;
    Corr(fInd,:) = out2.correlation;
    Ret1 = cat(1,out2.ret1,Ret1);
    Ret2 = cat(1,out2.ret2,Ret2);
    Ret3 = cat(1,out2.ret3,Ret3);
    Rel2 = cat(1,out2.rel2,Rel2);
    Rel3 = cat(1,out2.rel3,Rel3);
    if length(Ret1) ~= length(Ret3)
        break; % break if there is something wrong
    end
end
%%
mkdir(resultFolderName)
cd(resultFolderName)
save([resultFolderName,'DermingAndPearson'],'Slope','Corr');
save([resultFolderName,'Data'],'Ret1','Ret2','Ret3','Rel2','Rel3')
compoundPersonABplot(Ret1, Ret2, Ret3, Rel2, Rel3,0);
set(gcf,'Position',[200 200 1000,700])
saveas(gcf,[resultFolderName,outcomeName,'.bmp'])