function [MP1, MP2, cond] = mirrorPointExtract(S1, S2, pst)
% Compute histogram of the Stokes vectors to estimate the mirror state
%   mirrorPintExtract(S1, S2, pst)
%   will estmate the the mirror state (mirror state is equal to mirror
%   point)
%
%
%   S1,S2, the stokes vectors of the two channels. It can be bined or
%   non-bined.
%   pst, parameters for the assitance of mirror state estimation.
%
%   pst.makeOrth, determine whether make orthognal is needed, none by
%   default.  
%   pst.fwx, the size of kernel for stokes filtering.
%
%   pst.MP1, pst.MP2, the seed mirror state, which indicate the coarse area
%   of the true mirror state. The mirror state won't change largely for the
%   same sytem.
%
%   pst.searchRange determines the area of the search near the seed mirror
%   state. 10 by default. If no seed mirror state was provided, then a
%   large searchrange should be provided. such as 50.
%
%   MP1,MP2,  the estimated mirror state, the dimension depends on the input stokes vectors.
%
%   cond, the contrast of the mirror state. usually, the more prominent of
%   the estimated mirror state on the distribution map, the more reliable
%   is the mirror state estimation.
%
%   In this program, to speed up, a valid mirror state was provided.so the 
%   searchRange was set as 10. slightly change on the main program should 
%   be made if you want to generate a valid seed mirror state.
%
%
%   lastest updated at 9 July 2019,
%   by Qiaozhou Xiong, E150022@e.ntu.edu.sg, 
%   Nanyang Technological University,Singapore.
%
%%
m = size(S1);
if ~isfield(pst,'makeOrth')
    pst.makeOrth = 0;
end
h = filterGen(pst.fwx);
% h = 1;
if ~isfield(pst,'MP1')
    pst.MP1 = [0.6543 0.6967, 0.294];
    pst.MP2 = [-0.5287 0.8332, -0.1622];
end
if ~isfield(pst,'searchRange')
    pst.searchRange = 10;
end
searchRange = pst.searchRange;
% seed mirror state
MP1 = pst.MP1(m(3),:); MP2 = pst.MP2(m(3),:);
x1 = [asin(MP1(3)); atan2(MP1(2),MP1(1)).*cos(asin(MP1(3)))];
x2 = [asin(MP2(3)); atan2(MP2(2),MP2(1)).*cos(asin(MP2(3)))];

%% filtering
S1f = imfilter(S1,h,'circular');S2f = imfilter(S2,h,'circular');
I1 = imfilter(sqrt(sum(S1.^2,4)),h,'circular');I2 = imfilter(sqrt(sum(S2.^2,4)),h,'circular');
%% normalization
If1 = sqrt(dot(S1f,S1f,4));If2 = sqrt(dot(S2f,S2f,4));
dop1 = If1./I1;dop2 = If2./I2;
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
S1 = S1f./repmat(If1,[1,1,1,3]);
S2 = S2f./repmat(If2,[1,1,1,3]);
%% determine the threshold of the signal to use
mask = dop>.9;
mask(1:150,:) = 0;
mask(end-150:end,:) = 0;
if pst.makeOrth
    na = S1 + S2;
    nb = S1 - S2;
    S1 = na./repmat(sqrt(dot(na,na,4)),[1,1,1,3]);
    S2 = nb./repmat(sqrt(dot(nb,nb,4)),[1,1,1,3]);
end
%% Poincare shpere digitaliation
binsAz = linspace(-pi,pi,201);
binsEl = linspace(-pi/2,pi/2,101);
% determine the coordiate of the seed miror state
xEst1 = ceil((x1(1)-binsEl(1))/(binsEl(end)-binsEl(1))*numel(binsEl));
yEst1 = ceil((x1(2)-binsAz(1))/(binsAz(end)-binsAz(1))*numel(binsAz));

xEst2 = ceil((x2(1)-binsEl(1))/(binsEl(end)-binsEl(1))*numel(binsEl));
yEst2 = ceil((x2(2)-binsAz(1))/(binsAz(end)-binsAz(1))*numel(binsAz));

Sloc = S1;
elAngle = asin(Sloc(:,:,:,3));
azPos = atan2(Sloc(:,:,:,2),Sloc(:,:,:,1)).*cos(elAngle);
for wind = 1:m(3)
    elLoc = real(elAngle(:,:,wind));
    azLoc = real(azPos(:,:,wind));
    hh = hist2D(elLoc(mask),azLoc(mask),binsEl,binsAz);
    hh = imfilter(hh,filterGen(8)*(filterGen(10).'),'symmetric');
    HH1(:,:,wind) = hh;
end

Sloc = S2;
elAngle = asin(Sloc(:,:,:,3));
azPos = atan2(Sloc(:,:,:,2),Sloc(:,:,:,1)).*cos(elAngle);
for wind = 1:m(3)
    elLoc = real(elAngle(:,:,wind));
    azLoc = real(azPos(:,:,wind));
    hh = hist2D(elLoc(mask),azLoc(mask),binsEl,binsAz);
    hh = imfilter(hh,filterGen(8)*(filterGen(10).'),'symmetric');
    HH2(:,:,wind) = hh;
end

if m(3) == 9 %% determine whether is full-spectrum mirror state extraction
    for wind = 1:m(3)
        %determine the maximum density location of the stokes map at mirror
        %state cndidate of full-spectrum
        [maxDent,maxRow] = max(HH1(xEst1-searchRange:xEst1+searchRange,yEst1-searchRange:yEst1+searchRange,wind),[],1);
        [maxDent1,maxCol1] = max(maxDent,[],2);
        maxRow1 = maxRow(maxCol1);
        maxCol = maxCol1;
        maxRow = maxRow1;
        maxDent = maxDent1;
        yEstC1 = yEst1;
        xEstC1 = xEst1;
        cond(1,wind) = maxDent;
        maxColCh1(wind) = maxCol+yEstC1-searchRange-1;
        maxRowCh1(wind) = maxRow+xEstC1-searchRange-1;
        if maxDent<50
             maxColCh1(wind) = yEstC1;
             maxRowCh1(wind) = xEstC1;
        end
        %determine the maximum density location of the stokes map at mirror
        %state cndidate of full-spectrum
        [maxDent,maxRow] = max(HH2(xEst2-searchRange:xEst2+searchRange,yEst2-searchRange:yEst2+searchRange,wind),[],1);
        [maxDent2,maxCol2] = max(maxDent,[],2);
        maxRow2 = maxRow(maxCol2);
        maxCol = maxCol2;
        maxRow = maxRow2;
        maxDent = maxDent2;
        yEstC2 = yEst2;
        xEstC2 = xEst2;
        cond(2,wind) = maxDent;
        maxColCh2(wind) = maxCol+yEstC2-searchRange-1;
        maxRowCh2(wind) = maxRow+xEstC2-searchRange-1;
        if maxDent<50
             maxColCh2(wind) = yEstC2;
             maxRowCh2(wind) = xEstC2;
        end
    end
else %% full spectrum mirror state extraction
    wind = 1;
    %determine the maximum density location of the stokes map at mirror
    %state cndidate area 1
    [maxDent,maxRow] = max(HH1(xEst1-searchRange:xEst1+searchRange,yEst1-searchRange:yEst1+searchRange,wind),[],1);
    [maxDent1,maxCol1] = max(maxDent,[],2);
    maxRow1 = maxRow(maxCol1);
    %determine the maximum density location of the stokes map at mirror
    %state cndidate area 2
    [maxDent,maxRow] = max(HH1(xEst2-searchRange:xEst2+searchRange,yEst2-searchRange:yEst2+searchRange,wind),[],1);
    [maxDent2,maxCol2] = max(maxDent,[],2);
    maxRow2 = maxRow(maxCol2);
    maxCol = maxCol2;
    maxRow = maxRow2;
    maxDent = maxDent2;
    yEstC1 = yEst2;
    xEstC1 = xEst2;
    % which one to choose determine by the density
    if maxDent1>maxDent2
        maxCol = maxCol1;
        maxRow = maxRow1;
        maxDent = maxDent1;
        yEstC1 = yEst1;
        xEstC1 = xEst1;
    end
    cond(1,wind) = maxDent;
    maxColCh1(wind) = maxCol+yEstC1-searchRange-1;
    maxRowCh1(wind) = maxRow+xEstC1-searchRange-1;
    
    %determine the maximum density location of the stokes map at mirror
    %state cndidate area 1
    [maxDent,maxRow] = max(HH2(xEst1-searchRange:xEst1+searchRange,yEst1-searchRange:yEst1+searchRange,wind),[],1);
    [maxDent1,maxCol1] = max(maxDent,[],2);
    maxRow1 = maxRow(maxCol1);
    %determine the maximum density location of the stokes map at mirror
    %state cndidate area 2
    [maxDent,maxRow] = max(HH2(xEst2-searchRange:xEst2+searchRange,yEst2-searchRange:yEst2+searchRange,wind),[],1);
    [maxDent2,maxCol2] = max(maxDent,[],2);
    maxRow2 = maxRow(maxCol2);
    maxCol = maxCol2;
    maxRow = maxRow2;
    maxDent = maxDent2;
    yEstC2 = yEst2;
    xEstC2 = xEst2;
    % which one to choose determine by the density
    if maxDent1>maxDent2
        maxCol = maxCol1;
        maxRow = maxRow1;
        maxDent = maxDent1;
        yEstC2 = yEst1;
        xEstC2 = xEst1;
    end
    cond(2,wind) = maxDent;
    maxColCh2(wind) = maxCol+yEstC2-searchRange-1;
    maxRowCh2(wind) = maxRow+xEstC2-searchRange-1;
    
end
%
figure(2)
clf
imagesc(cat(1,HH1(:,:),HH2(:,:)))
hold on
plot(maxColCh1(:) + (0:numel(binsAz):(m(3)-1)*numel(binsAz))',maxRowCh1(:),'yo')
plot(maxColCh2(:) + (0:numel(binsAz):(m(3)-1)*numel(binsAz))',numel(binsEl) + maxRowCh2(:),'yo')
plot(yEstC1(:) + (0:numel(binsAz):(m(3)-1)*numel(binsAz))',xEstC1(:),'ro')
plot(yEstC2(:) + (0:numel(binsAz):(m(3)-1)*numel(binsAz))',numel(binsEl) + xEstC2(:),'ro')
title('Histogram of SOPs with estimated maxima')
xlabel('Azimuth for the different spectral bins')
ylabel('Elevation for S1 and S2, resp')
set(gca,'XTick',100:200:1800);set(gca,'XTicklabel',{'bin1','bin2','bin3','bin4','bin5','bin6','bin7','bin8','bin9'})
set(gca,'YTick',[1 90 110 200]);set(gca,'YTicklabel',{'-\pi/2','\pi/2','-\pi/2','\pi/2'})

%%
% use the maxima to define the MP
el = binsEl(maxRowCh1);
az = binsAz(maxColCh1);
MP1 = ([cos(el).*cos(az);cos(el).*sin(az);sin(el)]).';

el = binsEl(maxRowCh2);
az = binsAz(maxColCh2);
MP2 = ([cos(el).*cos(az);cos(el).*sin(az);sin(el)]).';