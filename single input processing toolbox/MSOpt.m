function out = MSOpt(S,x,pst,tauref,mask)
%w = mpSingleBin(S1,mp,dz) estimates the retardation vector given the 
%Stokes vectors S1 and the mirror point mp. dz is the axial range for
%filtering the resulting retardation vector.
dz = pst.dz;
Sf = imfilter(S,filterGen(pst.fwx).','circular');
If = sqrt(dot(Sf,Sf,4));
Sn = Sf./repmat(If,[1,1,1,3]);
sPlus = circshift(Sn,-dz);% Should be Nz by NAlines by 3
sMinus= circshift(Sn,dz);
MPnew = permute([cos(x(1)).*cos(x(2));cos(x(1)).*sin(x(2));sin(x(1))],[4,3,2,1]);

vp = sPlus-sMinus;
v = (sPlus+sMinus)/2;
dvu = MPnew-v;

dotVMP = v(:,:,:,1).*MPnew(:,:,:,1) + v(:,:,:,2).*MPnew(:,:,:,2) + v(:,:,:,3).*MPnew(:,:,:,3);
tau = bsxfun(@rdivide,cross(vp,dvu,4),1-dotVMP)/2/dz;
% ret = sqrt(dot(tau,tau,4));
% retref = sqrt(dot(tauref,tauref,4));
Rel = normalizeStokes3D(dvu);
mask = and(mask,(Rel>.5));
% corr = corrcoef(ret(mask),retref(mask));
% out = abs(corr(2,1)-1);
% out = sum(sum(bsxfun(@times,abs(ret-retref),mask),1),2);
% tau = imfilter(tau,filterGen(10)*(filterGen(10).'),'replicate');
% tauref = imfilter(tauref,filterGen(10)*(filterGen(10).'),'replicate');
out = sum(sum(bsxfun(@times,sum((tau-tauref).^2,4),mask),1),2);
% out = out./sum(sum(mask));