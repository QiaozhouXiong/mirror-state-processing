function out = PSProcessLocalSubsampling(S1,S2,procStruct)
%PSPROCESS - computes local birefringence of the Stokes vectors S1 and S2
% 
% The third dimension can either be the y coordinate, or the pc multiplexed
% dimension. This is indicated in the procStruct with the field pc.
% Using the dop is default for defining the correction mask
% Also, using the pc and pcEst are default
% Feb 2018: EuclideanRotation routine instead of SVD to speed up
% processing.

% the only two mandatory arguments
fwx = round(procStruct.fwx/2);
dz = procStruct.dz;
clipLimit = procStruct.clipLimit;

if ~isfield(procStruct,'fwz')
    fwz = 1;
else
    fwz = procStruct.fwz;
end

if isfield(procStruct,'fwy')
    fwy = procStruct.fwy;
else
    fwy = 1;
end

% check if pc processing along the third dimension
if isfield(procStruct,'pc')
    pc = procStruct.pc;
else
    pc = true;
end

% implementation of estimation of the relative rotation between the PAs of
% the subwindows
if isfield(procStruct,'pcEst')
    pcEst = procStruct.pcEst;
else
    pcEst = true; % pcEst is default
end

if isfield(procStruct,'dop')
    dopCheck = procStruct.dop;
else
    dopCheck = true; % dop is default
end

if isfield(procStruct,'dopTh')
    dopTh = procStruct.dopTh;
else
    dopTh = [0.6,1];
end

if isfield(procStruct,'dopu')
    dopuCheck = procStruct.dopu;
else
    dopuCheck = false;
end

if isfield(procStruct,'dopuTh')
    dopuTh = procStruct.dopuTh;
else
    dopuTh = [0.45,1];
end

% gaussian lateral filtering is default
if ~isfield(procStruct,'filter')
    procStruct.filter = 'gaussian';
end

if ~isfield(procStruct,'corrOut')
    corrOut = false;
else
    corrOut = procStruct.corrOut;
end

S1 = S1(:,1:2:end,:,:);
S2 = S2(:,1:2:end,:,:);

% padd Stokes vectors on both sides
prepad = floor(ceil(1.5*fwx)/2)+floor(fwx/2);
postpad = ceil(ceil(1.5*fwx)/2)+ceil(fwx/2);

S1 = cat(2,S1(:,end-prepad+1:end,:,:),S1,S1(:,1:postpad,:,:));
S2 = cat(2,S2(:,end-prepad+1:end,:,:),S2,S2(:,1:postpad,:,:));


I1 = sqrt(dot(S1,S1,4));
I2 = sqrt(dot(S2,S2,4));

if dopuCheck
    % pre-normalization for computation of DOPU
    Sn1 = S1./repmat(I1,[1,1,1,3]);
    Sn2 = S2./repmat(I2,[1,1,1,3]);
end

if ~pc && fwy>1
    % central slice if oop>1 and not pc processing
    i1 = squeeze(I1(:,:,ceil(fwy/2):end-ceil((fwy-1)/2)));
    i2 = squeeze(I2(:,:,ceil(fwy/2):end-ceil((fwy-1)/2)));
else
    i1 = I1;
    i2 = I2;
end

II = i1 + i2;

% filtering of S1 and S2
if isfield(procStruct,'filter') && strcmp(procStruct.filter,'mean') 
    S1 = imfilter(S1,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    S2 = imfilter(S2,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    I1 = imfilter(I1,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    I2 = imfilter(I2,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    if dopuCheck
        Sn1 = imfilter(Sn1,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
        Sn2 = imfilter(Sn2,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    end
elseif isfield(procStruct,'filter') && strcmp(procStruct.filter,'gaussian') 
    nx = (round(fwx*1.5)-1)/2;
    nz = (round(fwz*1.5)-1)/2;
    nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
    nz = linspace(-nz,nz,round(fwz*1.5))'*2*sqrt(log(2))/fwz;
    if fwz == 1
        nz = 0;
    end
    h = exp(-nz.^2)*exp(-nx.^2);
    h = h/sum(h(:));
    h = repmat(h,[1,1,fwy]);% out of plane it remains a simple mean filter
    S1 = imfilter(S1,h,'circular');
    S2 = imfilter(S2,h,'circular');
    I1 = imfilter(I1,h,'circular');
    I2 = imfilter(I2,h,'circular');
    if dopuCheck
        Sn1 = imfilter(Sn1,h,'circular');
        Sn2 = imfilter(Sn2,h,'circular');
    end
end


% final normalization
If1 = sqrt(dot(S1,S1,4));
If2 = sqrt(dot(S2,S2,4));

if dopCheck
    dop1 = If1./I1;
    dop2 = If2./I2;
    if pc
        dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);
    else
        dop = 1/2*(dop1+dop2);
    end
end

if dopuCheck
    dopu1 = sqrt(dot(Sn1,Sn1,4));
    dopu2 = sqrt(dot(Sn2,Sn2,4));
    if pc
        dopu = 1/2*mean(dopu1,3) + 1/2*mean(dopu2,3);
    else
        dopu = 1/2*(dopu1+dopu2);
    end
end

If = If1 + If2;

S1 = S1./repmat(If1,[1,1,1,3]);
S2 = S2./repmat(If2,[1,1,1,3]);

% force the two Stokes vectors to be orthogonal, which is equivalent to the
% lsq solution

% construct orthonormal tripod for these data points
na = S1 + S2;
nb = S1 - S2;
S1 = na./repmat(sqrt(dot(na,na,4)),[1,1,1,3]);
S2 = nb./repmat(sqrt(dot(nb,nb,4)),[1,1,1,3]);


% local birefringence analysis
S1plus = circshift(S1,-dz);
S1minus = circshift(S1,dz); 
S2plus = circshift(S2,-dz);
S2minus = circshift(S2,dz);   

% simple cross product of difference vectors to find PA
PA = cross(S1minus-S1plus,S2minus-S2plus,4);
PA = PA./repmat(max(sqrt(dot(PA,PA,4)),1e-9),[1,1,1,3]);

temp = dot(S1plus,PA,4).^2;% subtract tiny number to avoid division by zero in the next line
retSinW = real(acos((dot(S1plus,S1minus,4)-temp)./(1-temp + 1e-9)))/2/dz;
pm = sign((1-dot(S1plus,S1minus,4)).*(dot(S1plus-S1minus,S2plus+S2minus,4)));
PA = PA.*repmat(pm,[1,1,1,3]);
% 
if pc
    rpamean = mean(PA.*repmat(retSinW,[1,1,1,3]),3);
    rmean = sqrt(dot(rpamean,rpamean,4));
    out.stdpa = 1-rmean(:,prepad+1:end-postpad)./mean(abs(retSinW(:,prepad+1:end-postpad,:)),3);    
    out.rmean = rmean(:,prepad+1:end-postpad)*100/4.8*180/pi;
    %retSinW = mean(abs(retSinW),3);
end

% temp
Wcorr = PA.*repmat(retSinW,[1,1,1,3]);
out.tauref = Wcorr(:,prepad+1:end-postpad,:,:);
dim = size(PA);
mask = (dop>.8);

mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;

mid = ceil(size(retSinW,3)/2);
B = bsxfun(@times,Wcorr(:,:,mid,:),mask);

for wind = [(1:mid-1),(mid+1:dim(3))]
    A = bsxfun(@times,Wcorr(:,:,wind,:),mask);

    C(1,1) = sum(sum(A(:,:,1,1).*B(:,:,1,1)));
    C(2,1) = sum(sum(A(:,:,1,1).*B(:,:,1,2)));
    C(3,1) = sum(sum(A(:,:,1,1).*B(:,:,1,3)));

    C(1,2) = sum(sum(A(:,:,1,2).*B(:,:,1,1)));
    C(2,2) = sum(sum(A(:,:,1,2).*B(:,:,1,2)));
    C(3,2) = sum(sum(A(:,:,1,2).*B(:,:,1,3)));

    C(1,3) = sum(sum(A(:,:,1,3).*B(:,:,1,1)));
    C(2,3) = sum(sum(A(:,:,1,3).*B(:,:,1,2)));
    C(3,3) = sum(sum(A(:,:,1,3).*B(:,:,1,3)));

    R = reshape(euclideanRotation(C(:)),[3,3]);

    temp1 = R(1,1)*Wcorr(:,:,wind,1) + R(1,2)*Wcorr(:,:,wind,2) + R(1,3)*Wcorr(:,:,wind,3);
    temp2 = R(2,1)*Wcorr(:,:,wind,1) + R(2,2)*Wcorr(:,:,wind,2) + R(2,3)*Wcorr(:,:,wind,3);
    temp3 = R(3,1)*Wcorr(:,:,wind,1) + R(3,2)*Wcorr(:,:,wind,2) + R(3,3)*Wcorr(:,:,wind,3);

    Wcorr(:,:,wind,1) = temp1;
    Wcorr(:,:,wind,2) = temp2;
    Wcorr(:,:,wind,3) = temp3;
    rcorrglobal(:,wind) = decomposeRot(R(:));
end

rpamean = mean(Wcorr,3);
rmeancorr2 = sqrt(dot(rpamean,rpamean,4))*0.0239;


if pcEst
    % implementation of estimation of relative rotations between the
    % different subwindows

    PA = PA.*repmat(retSinW,[1,1,1,3]);
    Wcorr = PA;
    dim = size(PA);

    if dopCheck
        mask = (dop>dopTh(1)).*(dop<=dopTh(2));
    elseif dopuCheck
        mask = (dopu>dopuTh(1)).*(dopu<=dopuTh(2));
    end

    mid = ceil(size(retSinW,3)/2);
    Wcorr = Wcorr.*repmat(mask,[1,1,dim(3),3]);

    ref = Wcorr(:,:,mid,:);
    C = zeros([3,3,dim(2)]);

    h = ones(1,fwx)/3;
    rcorrconv = zeros(3,size(C,3),size(retSinW,3));
    for wind = [(1:mid-1),(mid+1:dim(3))]
        C(1,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,1),1),h,'same');
        C(2,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,2),1),h,'same');
        C(3,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,3),1),h,'same');

        C(1,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,1),1),h,'same');
        C(2,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,2),1),h,'same');
        C(3,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,3),1),h,'same');

        C(1,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,1),1),h,'same');
        C(2,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,2),1),h,'same');
        C(3,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,3),1),h,'same');

        R = reshape(euclideanRotation(reshape(C,[9,size(C,3)])),[3,3,size(C,3)]);

        temp1 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(1,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(1,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(1,3,:),1),[dim(1),1]);
        temp2 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(2,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(2,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(2,3,:),1),[dim(1),1]);
        temp3 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(3,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(3,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(3,3,:),1),[dim(1),1]);

        Wcorr(:,:,wind,1) = temp1;
        Wcorr(:,:,wind,2) = temp2;
        Wcorr(:,:,wind,3) = temp3;
        
        rcorrconv(:,:,wind) = decomposeRot(reshape(R,[9,size(C,3)]));
    end
    
    rpamean = mean(Wcorr,3);
    rmeancorr = sqrt(dot(rpamean,rpamean,4))*0.0239;

end

if pcEst
    out.rmeancorr = rmeancorr(:,prepad+1:end-postpad);
    out.PAcorr = rpamean(:,prepad+1:end-postpad,:,:)./repmat(rmeancorr(:,prepad+1:end-postpad),[1,1,1,3]);
    %    out.errCorr = errCorr;
end

out.I = squeeze(II(:,prepad+1:end-postpad,:));
out.If = squeeze(mean(I1(:,prepad+1:end-postpad,:)+I1(:,prepad+1:end-postpad,:),3));
out.ret = squeeze(retSinW(:,prepad+1:end-postpad,:))*100/4.8*180/pi;
out.PA = squeeze(PA(:,prepad+1:end-postpad,:,:));
if dopuCheck
    out.dopu1 = dopu1(:,prepad+1:end-postpad,:);
    out.dopu2 = dopu2(:,prepad+1:end-postpad,:);
    out.dopu = dopu(:,prepad+1:end-postpad,:);
end
    
if dopCheck
    out.dop1 = dop1(:,prepad+1:end-postpad,:);
    out.dop2 = dop2(:,prepad+1:end-postpad,:);
    out.dop = dop(:,prepad+1:end-postpad,:);
end

out.rcorrconv = rcorrconv(:,prepad+1:end-postpad,:);
out.rcorrglobal = rcorrglobal;
out.rmeancorr2 = rmeancorr2(:,prepad+1:end-postpad);


