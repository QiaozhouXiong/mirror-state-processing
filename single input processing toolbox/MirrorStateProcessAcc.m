function out = MirrorStateProcessAcc(S1,S2,pst)
% Compute local retardation with single channel and generate a compound
% estimated reliability.
%
% S1,S2: the Stokes vectors of the two polarization channels
%
% pst.MP1, pst.MP2, mirror state of the two polarization channels
%
% pst.fwx, kernel size for lateral fitering
%
% pst.dz, the axial shift applied to stokes vector for rotation vector
% reasoning.
%
% pst.clipLimit, the lower surface of the sheath, removed from the analysis
%
%
% laste updated at 10 July 2o19 by Qiaozhou Xiong (E150022@e.ntu.edu.sg)
%% mirror point prepartion
dim = size(S1);
MP1(1,1,:,:) = pst.MP1;
MP1 = repmat(MP1,[1024, dim(2), 1, 1]);
MP2(1,1,:,:) = pst.MP2;
MP2 = repmat(MP2,[1024, dim(2), 1, 1]);
%%
%-define a threshold to mask the area for the correction matrix calculation
MPAthre = 0.5;
fwx = pst.fwx;
dz = pst.dz;
clipLimit = pst.clipLimit;
h = filterGen(fwx).';
S1f = imfilter(S1,h,'circular');
S2f = imfilter(S2,h,'circular');
I1 = imfilter(sqrt(sum(S1.^2,4)),h,'circular');
I2 = imfilter(sqrt(sum(S2.^2,4)),h,'circular');

If1 = sqrt(dot(S1f,S1f,4));
If2 = sqrt(dot(S2f,S2f,4));
%% estimate the SNR
SNRmetric1 = bsxfun(@minus,10*log10(I1),mean(mean(10*log10(I1(600:800,:,:)),1),2));
SNRmetric2 = bsxfun(@minus,10*log10(I2),mean(mean(10*log10(I2(600:800,:,:)),1),2));

% dop calculation
dop1 = If1./I1;
dop2 = If2./I2;
% normalization
S1n = S1f./repmat(If1,[1,1,1,3]);
S2n = S2f./repmat(If2,[1,1,1,3]);
% Generate mask, based on DOP, and excluding the first and last depth pixels.
mask = mean(dop1,3)>.8;
mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;
diat = dot(squeeze(S1n(:,:,5,:)),squeeze(S2n(:,:,5,:)),3);
diat(~mask) = 0;
out.diat = diat; % for the diattenuation check
%% mirror point reconstruction 1 we used more accurate method here
sPlus = circshift(S1n,-dz);% Should be Nz by NAlines by 3
sMinus= circshift(S1n,dz);
vp = sPlus-sMinus;
dvu = MP1-sPlus;
tau = cross(vp,dvu,4);
%% normalize 
tau = tau./sqrt(sum(tau.^2,4));
ret = atan2(sum(cross(sPlus,sMinus,4).*tau,4),(sum(sPlus.*sMinus,4)-sum(tau.*sPlus,4).^2));
tau = tau.*ret;
sina = sqrt(sum((cross(tau,sPlus,4)./sqrt(sum(tau.^2,4))./sqrt(sum(sPlus.^2,4))).^2,4));
%% There are three ways to qunatify the distance
relW = (normalizeStokes(sPlus-MP1) + normalizeStokes(sMinus-MP1))/2;
% relW = normalizeStokes(dvu);
% select area with good reliability to calcualte rotation matrix
maskMpa = relW>MPAthre;
% now, align the spectral bins
% mask out regions with low reliability
dim = size(tau);
mid = ceil(dim(3)/2);
maskDop = dop1>0.8;
out.tau = tau;
out.relW1 = relW;
out.SNR1 = SNRmetric1;
out.SNR2 = SNRmetric2;
out.binRet1 = normalizeStokes(tau);
taucorr = tau;
if strcmp(pst.relm, 'dist')
%     SNRmetric1(SNRmetric1<0) = 0;
%     SNRmetric1(SNRmetric1>25) = 25;
    avgWeg = (polyval(pst.Apval,SNRmetric1).*relW+polyval(pst.Bpval,SNRmetric1)).^2;
%     avgWeg = avgWeg.*(relW>0.2);
else
    avgWeg = (polyval(pst.Apval,SNRmetric1).*sina+polyval(pst.Bpval,SNRmetric1)).^2;
end
%% weight for spectral alignment
weight = avgWeg.*maskDop.*maskMpa.*(SNRmetric1>10).*mask;
% weight = avgWeg;
out.avgWeg1 = avgWeg;
B = squeeze(taucorr(:,:,mid,:));
for wind = [(1:mid-1),(mid+1:dim(3))]
    ww = weight(:,:,wind).*weight(:,:,mid);
    A = squeeze(bsxfun(@times,taucorr(:,:,wind,:),ww));
    
    C(1,1) = sum(sum(A(:,:,1).*B(:,:,1)));
    C(2,1) = sum(sum(A(:,:,1).*B(:,:,2)));
    C(3,1) = sum(sum(A(:,:,1).*B(:,:,3)));
    
    C(1,2) = sum(sum(A(:,:,2).*B(:,:,1)));
    C(2,2) = sum(sum(A(:,:,2).*B(:,:,2)));
    C(3,2) = sum(sum(A(:,:,2).*B(:,:,3)));
    
    C(1,3) = sum(sum(A(:,:,3).*B(:,:,1)));
    C(2,3) = sum(sum(A(:,:,3).*B(:,:,2)));
    C(3,3) = sum(sum(A(:,:,3).*B(:,:,3)));
    
    R = reshape(euclideanRotation(C(:)),[3,3]);
    
    temp1 = R(1,1)*taucorr(:,:,wind,1) + R(1,2)*taucorr(:,:,wind,2) + R(1,3)*taucorr(:,:,wind,3);
    temp2 = R(2,1)*taucorr(:,:,wind,1) + R(2,2)*taucorr(:,:,wind,2) + R(2,3)*taucorr(:,:,wind,3);
    temp3 = R(3,1)*taucorr(:,:,wind,1) + R(3,2)*taucorr(:,:,wind,2) + R(3,3)*taucorr(:,:,wind,3);
    
    taucorr(:,:,wind,1) = temp1;
    taucorr(:,:,wind,2) = temp2;
    taucorr(:,:,wind,3) = temp3;
    rcorrglobal(:,wind) = decomposeRot(R(:));
end
avgWeg(avgWeg<=0) = 1e-9;
taumean = squeeze(sum(bsxfun(@times,taucorr,avgWeg),3));
taumean = bsxfun(@rdivide, taumean,sum(avgWeg,3));
taumean(isnan(taumean)) = 0;
tauErr = sqrt(sum((bsxfun(@minus,taucorr,permute(taumean,[1 2 4 3]))).^2,4));
tauErr(repmat(mask,[1 1 3])) = 0;
out.binError1 = sum(squeeze(sum(tauErr)));

out.rcorrglobal1 = rcorrglobal;
ret1 = sqrt(sum(taumean.^2,3))*0.0239/2/dz;
out.ret1 = ret1;
%% reconstruction with channel 2
sPlus = circshift(S2n,-dz);% Should be Nz by NAlines by 3
sMinus= circshift(S2n,dz);
vp = sPlus-sMinus;

dvu = MP2-sPlus;
tau = cross(vp,dvu,4);
taunorm = sqrt(sum(tau.^2,4));
tau = tau./taunorm;
ret = atan2(sum(cross(sPlus,sMinus,4).*tau,4),(sum(sPlus.*sMinus,4)-sum(tau.*sPlus,4).^2));
tau = tau.*ret;
%%
mask = mean(dop2,3)>.8;
mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;
sina = sqrt(sum((cross(tau,sPlus,4)./sqrt(sum(tau.^2,4))./sqrt(sum(sPlus.^2,4))).^2,4));
% distance as the mean of the distance to mirror state respectively
relW = (normalizeStokes(sPlus-MP2) + normalizeStokes(sMinus-MP2))/2;
maskMpa = relW>MPAthre;
% mask out regions with low reliability
dim = size(tau);
mid = ceil(dim(3)/2);
maskDop = dop2>0.8;
taucorr = tau;
out.binRet2 = normalizeStokes(tau);
if strcmp(pst.relm, 'dist')
    SNRmetric2(SNRmetric2<0) = 0;
    SNRmetric2(SNRmetric2>25) = 25;
    avgWeg = (polyval(pst.Apval,SNRmetric2).*relW+polyval(pst.Bpval,SNRmetric2)).^2;
%     avgWeg = avgWeg.*(relW>0.2);
else
    avgWeg = (polyval(pst.Apval,SNRmetric2).*sina+polyval(pst.Bpval,SNRmetric2)).^2;
end
weight = avgWeg.*maskDop.*maskMpa.*(SNRmetric2>10).*mask;
B = squeeze(taucorr(:,:,mid,:));
for wind = [(1:mid-1),(mid+1:dim(3))]
    ww = weight(:,:,wind).*weight(:,:,mid);
    A = squeeze(bsxfun(@times,taucorr(:,:,wind,:),ww));
    
    C(1,1) = sum(sum(A(:,:,1).*B(:,:,1)));
    C(2,1) = sum(sum(A(:,:,1).*B(:,:,2)));
    C(3,1) = sum(sum(A(:,:,1).*B(:,:,3)));
    
    C(1,2) = sum(sum(A(:,:,2).*B(:,:,1)));
    C(2,2) = sum(sum(A(:,:,2).*B(:,:,2)));
    C(3,2) = sum(sum(A(:,:,2).*B(:,:,3)));
    
    C(1,3) = sum(sum(A(:,:,3).*B(:,:,1)));
    C(2,3) = sum(sum(A(:,:,3).*B(:,:,2)));
    C(3,3) = sum(sum(A(:,:,3).*B(:,:,3)));
    
    R = reshape(euclideanRotation(C(:)),[3,3]);
    
    temp1 = R(1,1)*taucorr(:,:,wind,1) + R(1,2)*taucorr(:,:,wind,2) + R(1,3)*taucorr(:,:,wind,3);
    temp2 = R(2,1)*taucorr(:,:,wind,1) + R(2,2)*taucorr(:,:,wind,2) + R(2,3)*taucorr(:,:,wind,3);
    temp3 = R(3,1)*taucorr(:,:,wind,1) + R(3,2)*taucorr(:,:,wind,2) + R(3,3)*taucorr(:,:,wind,3);
    
    taucorr(:,:,wind,1) = temp1;
    taucorr(:,:,wind,2) = temp2;
    taucorr(:,:,wind,3) = temp3;
    rcorrglobal(:,wind) = decomposeRot(R(:));
end
avgWeg(avgWeg<=0) = 1e-9;
taumean = squeeze(sum(bsxfun(@times,taucorr,avgWeg),3));
taumean = bsxfun(@rdivide, taumean,sum(avgWeg,3));
taumeant = squeeze(mean(taucorr,3));
nanMask = isnan(taumean);
taumean(nanMask) = taumeant(nanMask);
tauErr = sqrt(sum((bsxfun(@minus,taucorr,permute(taumean,[1 2 4 3]))).^2,4));
tauErr(repmat(mask,[1 1 3])) = 0;
out.binError2 = sum(squeeze(sum(tauErr)));

out.rcorrglobal2 = rcorrglobal;
out.avgWeg2 = avgWeg;
out.relW2 = relW;
ret2 = sqrt(sum(taumean.^2,3))*0.0239/2/dz;
% ret2(sum(avgWeg,3)<30) = 0;
out.ret2 = ret2;
%% differenciate which channle
if dot(pst.MP1(1,:),[1 0 0])<dot(pst.MP2(1,:),[1 0 0])
    temp1 = out.binError2;
    temp2 = out.rcorrglobal2;
    temp3 = out.avgWeg2;
    out.binError2 = out.binError1;
    out.rcorrglobal2 = out.rcorrglobal1;
    out.avgWeg2 = out.avgWeg1;
    out.ret2 = ret1;
    
    out.binError1 = temp1;
    out.rcorrglobal1 = temp2;
    out.avgWeg1 = temp3;
    out.ret1 = ret2;
end