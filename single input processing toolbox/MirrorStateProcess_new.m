function out = MirrorStateProcess_new(S1,S2,pst)
if ~isfield(pst,'makeOrth')
    pst.makeOrth = 0;
end
dim = size(S1);
MP1(1,1,:,:) = pst.MP1;
MP1 = repmat(MP1,[1024, dim(2), 1, 1]);
MP2(1,1,:,:) = pst.MP2;
MP2 = repmat(MP2,[1024, dim(2), 1, 1]);
MPAthre = 0.4;
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

SNRmetric1 = 10*log10(I1)-60;
SNRmetric2 = 10*log10(I2)-60;
SNRmetric1(SNRmetric1>20) = 20;
SNRmetric2(SNRmetric2>20) = 20;
SNRmetric1(SNRmetric1<0) = 0;
SNRmetric2(SNRmetric2<0) = 0;
SNRmetric1 = 10.^(SNRmetric1./20);
SNRmetric2 = 10.^(SNRmetric2./20);

dop1 = If1./I1;
dop2 = If2./I2;
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
S1n = S1f./repmat(If1,[1,1,1,3]);
S2n = S2f./repmat(If2,[1,1,1,3]);
if pst.makeOrth
    na = S1n + S2n;
    nb = S1n - S2n;
    S1n = na./repmat(sqrt(dot(na,na,4)),[1,1,1,3]);
    S2n = nb./repmat(sqrt(dot(nb,nb,4)),[1,1,1,3]);
end
% Generate mask, based on DOP, and excluding the first and last depth
% pixels.
mask = dop>.8;
mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;
diat = dot(squeeze(S1n(:,:,5,:)),squeeze(S2n(:,:,5,:)),3);
diat(~mask) = 0;
out.diat = diat;
%% mirror point reconstruction 1
sPlus = circshift(S1n,-dz);% Should be Nz by NAlines by 3
sMinus= circshift(S1n,dz);

vp = sPlus-sMinus;
v = (sPlus+sMinus)/2;
dvu = MP1-v;

dotVMP = v(:,:,:,1).*MP1(:,:,:,1) + v(:,:,:,2).*MP1(:,:,:,2) + v(:,:,:,3).*MP1(:,:,:,3);
tau = bsxfun(@rdivide,cross(vp,dvu,4),1-dotVMP)/2/dz;
sina = sqrt(sum((cross(tau,v,4)./sqrt(sum(tau.^2,4))./sqrt(sum(v.^2,4))).^2,4));
relW = (sqrt(sum(dvu.^2,4)));%1-dotVMP.^2;
maskMpa = sina>MPAthre;
% tau should match tauref

% now, align the spectral bins
% mask out regions with low reliability
dim = size(tau);
mid = ceil(dim(3)/2);
maskDop = dop1>0.8;
out.tau = tau;
taucorr = tau;
load('E:\Single input intravascular polarimetry with mirror state\Step 3 Single input processing\relMet.mat');
avgWeg = (polyval(Apval,SNRmetric1).*sina+polyval(Bpval,SNRmetric1)).^2;
weight = avgWeg.*maskDop.*maskMpa.*(SNRmetric1>2).*mask;
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
ret1 = sqrt(sum(taumean.^2,3));
out.ret1 = ret1*0.0239;
%% reconstruction with channel 2
sPlus = circshift(S2n,-dz);% Should be Nz by NAlines by 3
sMinus= circshift(S2n,dz);

vp = sPlus-sMinus;
v = (sPlus+sMinus)/2;
dvu = MP2-v;

dotVMP = v(:,:,:,1).*MP2(:,:,:,1) + v(:,:,:,2).*MP2(:,:,:,2) + v(:,:,:,3).*MP2(:,:,:,3);
tau = bsxfun(@rdivide,cross(vp,dvu,4),1-dotVMP)/2/dz;
sina = sqrt(sum((cross(tau,v,4)./sqrt(sum(tau.^2,4))./sqrt(sum(v.^2,4))).^2,4));
relW = (sqrt(sum(dvu.^2,4)));%1-dotVMP.^2;
maskMpa = sina>MPAthre;
% tau should match tauref

% now, align the spectral bins
% mask out regions with low reliability
dim = size(tau);
mid = ceil(dim(3)/2);
maskDop = dop2>0.8;
taucorr = tau;
% avgWeg = ((1.3e-4*SNRmetric2.^3-0.0716*SNRmetric2.^2+1.61*SNRmetric2-1.943).*sina-6e-4*SNRmetric2.^2+0.00919*SNRmetric2+11.8645).^2;
% relW = (sqrt(sum(dvu.^2,4)));%1-dotVMP.^2;
% avgWeg = exp(-(0.21*SNRmetric2-2.47).*(-0.49*relW+0.195)).*(0.815*SNRmetric2-3.5).*relW;
avgWeg = ((polyval(Apval,SNRmetric2).*sina+polyval(Bpval,SNRmetric2))).^2;
avgWeg(avgWeg<0) = 0;
weight = avgWeg.*maskDop.*maskMpa.*(SNRmetric2>2).*mask;
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
ret2 = sqrt(sum(taumean.^2,3));
% ret2(sum(avgWeg,3)<30) = 0;
out.ret2 = ret2*0.0239;
%% differenciate which channle
if dot(pst.MP1(1,:),[1 0 0])<dot(pst.MP2(1,:),[1 0 0])
    temp1 = out.binError2;
    temp2 = out.rcorrglobal2;
    temp3 = out.avgWeg2;
    out.binError2 = out.binError1;
    out.rcorrglobal2 = out.rcorrglobal1;
    out.avgWeg2 = out.avgWeg1;
    out.ret2 = ret1*0.0239;
    
    out.binError1 = temp1;
    out.rcorrglobal1 = temp2;
    out.avgWeg1 = temp3;
    out.ret1 = ret2*0.0239;
end