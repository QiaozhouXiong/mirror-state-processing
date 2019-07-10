
% Script to compute the PDFs of the measured retardation vectors as a
% function of SNR and position relative to the mirror point
clear;
pathdef;
Nrep = 10000; % number of repetitions per SNR
SNR = linspace(6,25,41); % log scale
% SNRL = SNR;
SNRL = power(10, SNR/20);% linear scale
Nens = 10; % number of points averaged
Nnoise = 40; % number of repeats for each noise realization
% dret = linspace(0,0.8,Nrep); % range of retardation
% rand('state',0); 
dret = sqrt(2*log2(1./rand(1,Nrep)))*0.24;
% dret = 0.8;
mode = 'accurate';
% mode = 'approx';

clear WW WWGT Sina Theta;
for snrInd = 1:numel(SNR)
    snr = SNRL(snrInd);

    beta = randn(3,Nrep);
    betadiff = beta./sqrt(sum(beta.^2,1)).*dret; % rotation vector between two meausured Stokes vectors
    betamp = beta./sqrt(sum(beta.^2,1)).*rand(1,Nrep)*2*pi; % rotation vector from first Stokes vector to mirror point

    s1 = randn(3,Nrep);
    s1 = s1./sqrt(sum(s1.^2,1));

    RR = makeRot(betadiff);
    s2 = RR(1:3,:).*s1(1,:) + RR(4:6,:).*s1(2,:) + RR(7:9,:).*s1(3,:);

    RR = makeRot(betamp);
    mp = RR(1:3,:).*s1(1,:) + RR(4:6,:).*s1(2,:) + RR(7:9,:).*s1(3,:);

    % convert Stokes vectors back to Jones vectors for addition of noise
    phi = atan2(s1(3,:),s1(2,:));
    theta = atan2(sqrt(s1(2,:).^2 + s1(3,:).^2),s1(1,:));
    e1 = cat(1,cos(theta/2),sin(theta/2).*exp(1i*phi));
    
    phi = atan2(s2(3,:),s2(2,:));
    theta = atan2(sqrt(s2(2,:).^2 + s2(3,:).^2),s2(1,:));
    e2 = cat(1,cos(theta/2),sin(theta/2).*exp(1i*phi));

    % Speckle amplitude
    amp1 = abs(randn(1,Nrep,Nnoise,Nens) + 1i*randn(1,Nrep,Nnoise,Nens));
    amp2 = abs(randn(1,Nrep,Nnoise,Nens) + 1i*randn(1,Nrep,Nnoise,Nens));

    % noise amplitude
    noise1 = randn(2,Nrep,Nnoise,Nens) + 1i*randn(2,Nrep,Nnoise,Nens);
    noise2 = randn(2,Nrep,Nnoise,Nens) + 1i*randn(2,Nrep,Nnoise,Nens);
    
    eloc = e1.*amp1 + noise1/snr;
    S1 = mean(cat(1,abs(eloc(1,:,:,:)).^2-abs(eloc(2,:,:,:)).^2,2*real(eloc(1,:,:,:).*conj(eloc(2,:,:,:))),-2*imag(eloc(1,:,:,:).*conj(eloc(2,:,:,:)))),4);
    
    eloc = e2.*amp2 + noise2/snr;
    S2 = mean(cat(1,abs(eloc(1,:,:,:)).^2-abs(eloc(2,:,:,:)).^2,2*real(eloc(1,:,:,:).*conj(eloc(2,:,:,:))),-2*imag(eloc(1,:,:,:).*conj(eloc(2,:,:,:)))),4);

    S1 = S1./sqrt(sum(S1.^2,1));
    S2 = S2./sqrt(sum(S2.^2,1));

    if strcmp(mode,'accurate')
        % Accurate formulation
        vp = S2-S1;
        v = (S1+S2)/2;
        dvu = mp-v;

        %tau = cross(dvu,vp,1);
        tau = cross(vp,dvu,1);
        taunorm = sqrt(sum(tau.^2,1));
        tau = tau./taunorm;
        ret = atan2(sum(cross(S1,S2,1).*tau,1),(sum(S1.*S2,1)-sum(tau.*S1,1).^2));
        W = tau.*ret;
    elseif strcmp(mode,'approx')
        vp = S2-S1;
        v = (S1+S2)/2;
        dvu = mp-v;

        W = cross(vp,dvu,1)./(sum(v.^2,1)-sum(v.*mp,1));
        ret = sqrt(sum(W.^2,1));
    end
    
    WWGT(:,:,1,snrInd) = betadiff;
    WW(:,:,:,snrInd) = W;
    D(:,:,snrInd) = squeeze(sqrt(sum(dvu.^2,1)));
    D1(:,:,snrInd) = squeeze(sqrt(sum((mp-S1).^2,1)));
    D2(:,:,snrInd) = squeeze(sqrt(sum((mp-S2).^2,1)));
    Sina(:,:,snrInd) = sqrt(sum((cross((S1+S2)/2,W)./sqrt(sum(((S1+S2)/2).^2,1))./sqrt(sum(W.^2,1))).^2,1));
    

snrInd
end
WWGT = repmat(WWGT,[1,1,Nnoise,1]);
errRetV = squeeze(sqrt(sum((WW-WWGT).^2,1)));
errRet = squeeze(sqrt(sum(WW.^2,1))-sqrt(sum(WWGT.^2,1)));
%dim 1: changing ret and D
%dim 2: noise repeat
%dim 3: SNR

errV = WW-WWGT;

RetGT = squeeze(sqrt(sum(WWGT.^2,1)));
Ret = squeeze(sqrt(sum(WW.^2,1)));

%%
Dbins = linspace(0.01,1.8,81);
errSbins = linspace(-1,1,101);

for snrInd = 1:numel(SNR)
    
    errWloc = errRet(:,:,snrInd);
    
    errQloc = squeeze(errV(1,:,:,snrInd));
    errUloc = squeeze(errV(2,:,:,snrInd));
    errVloc = squeeze(errV(3,:,:,snrInd));
    
    
%     Dloc = min(D1(:,:,snrInd), D2(:,:,snrInd));

%     Dloc = sqrt(D1(:,:,snrInd).*D2(:,:,snrInd));
    Dloc = D(:,:,snrInd);
    [dsorted,inds] = sort(Dloc(:),'ascend');
    
    indStart = 1;
    for indd = 1:numel(Dbins)-1
        indEnd = find(dsorted<Dbins(indd+1),1,'last');
           
        NN(indd,snrInd) = indEnd-indStart+1;
        temp = errQloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHQ(:,indd,snrInd) = hh/sum(hh);
        meanQErr(indd,snrInd) = mean(temp);
        sigQErr(indd,snrInd) = std(temp);

        temp = errUloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHU(:,indd,snrInd) = hh/sum(hh);
        meanUErr(indd,snrInd) = mean(temp);
        sigUErr(indd,snrInd) = std(temp);

        temp = errVloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHV(:,indd,snrInd) = hh/sum(hh);
        meanVErr(indd,snrInd) = mean(temp);
        sigVErr(indd,snrInd) = std(temp);

        indStart = indEnd+1;
    end
    snrInd
end
HH = (HHQ + HHU + HHV)/3;
meanErr = (meanQErr + meanUErr + meanVErr)/3;
sigErr = (sigQErr + sigUErr + sigVErr)/3;
%%

figure(11)
clf
snrInd = 10;
subplot(2,3,1)
imagesc(Dbins(1:end-1),errSbins,HH(:,:,snrInd))
xlabel('D')
ylabel('errRet')
title(sprintf('errQUV, SNR = %.1f',SNR(snrInd)))
colorbar


subplot(2,3,2)
imagesc(SNR,Dbins(1:end-1),meanErr,[-.1,.1])
xlabel('SNR')
ylabel('D')
colorbar
title(sprintf('meanErrQUV'))


subplot(2,3,3)
imagesc(SNR,Dbins(1:end-1),sigErr)
xlabel('SNR')
ylabel('D')
colorbar
title(sprintf('stdErrQUV'))


dInd = 40;
subplot(2,3,4)
plot(errSbins,HHQ(:,dInd,snrInd))
hold on
plot(errSbins,HHU(:,dInd,snrInd))
plot(errSbins,HHV(:,dInd,snrInd))
xlabel('Error')
title(sprintf('At SNR = %.2f and D = %.2f-%.2f',SNR(snrInd),Dbins(dInd),Dbins(dInd+1)))
legend('Q','U','V')


subplot(2,3,5)
plot(Dbins(1:end-1),meanErr')
xlabel('D')
ylabel('Ret')
title(sprintf('Mean error'))

subplot(2,3,6)
plot(Dbins(1:end-1),sigErr')
xlabel('D')
ylabel('Ret')
title(sprintf('Std error'))

%%
%%
figure(4);clf; snrInd = 20;
subplot(2,2,1)
imagesc(Dbins(1:end-1),errSbins,HH(:,:,snrInd))
xlabel('D')
ylabel('errRet');axis xy;
title(sprintf('errQUV, SNR = %.1f',SNR(snrInd)))
colorbar  

subplot(2,2,2)
imagesc(SNR,Dbins(10:end-1),sigErr(10:end,:))
xlabel('SNR')
ylabel('D')
colorbar; axis xy;
title(sprintf('stdErrQUV'))

dInd = 40;
subplot(2,2,3)
plot(errSbins,HHQ(:,snrInd,dInd),'linewidth',2)
hold on
plot(errSbins,HHU(:,snrInd,dInd),'linewidth',2)
plot(errSbins,HHV(:,snrInd,dInd),'linewidth',2)
xlabel('Error');ylabel('Counts')
set(gca,'linewidth',1)
title(sprintf('At SNR = %.2f and D = %.2f-%.2f', SNR(snrInd),Dbins(dInd),Dbins(dInd+1)))
legend('Q','U','V');legend BOXOFF

subplot(2,2,4)
h1 = plot(Dbins(12:end-1),sigQErr(12:end, 10:9:30),'linewidth',2);
xlabel('D');set(gca,'linewidth',1)
ylabel('Standard deviation')
legend(['SNR=', num2str(SNR(10),2),'dB'],...
    ['SNR=', num2str(SNR(19),2),'dB'],...
    ['SNR=', num2str(SNR(28),2),'dB']);
legend BOXOFF
title(sprintf('Std error'))
set(gcf,'color','w');set(gca,'color','w')
axis([0 2, 0, 4]);
%% modeling
distanceModeling
axis([0 2, 0, 4]);
% Additional analysis in terms of Sina and Theta, instead of D
%%
Sinabins = linspace(0,1,81);
errSbins = linspace(-1,1,101);

for snrInd = 1:numel(SNR)
    
    errQloc = squeeze(errV(1,:,:,snrInd));
    errUloc = squeeze(errV(2,:,:,snrInd));
    errVloc = squeeze(errV(3,:,:,snrInd));

    Dloc = (Sina(:,:,snrInd));
    [dsorted,inds] = sort(Dloc(:),'ascend');
    
    indStart = 1;
    for indd = 1:numel(Sinabins)-1
        indEnd = find(dsorted<Sinabins(indd+1),1,'last');
        
        temp = Dloc(inds(indStart:indEnd));
        
        NN(indd,snrInd) = indEnd-indStart+1;
        temp = errQloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHQ(:,indd,snrInd) = hh/sum(hh);
        meanQErr(indd,snrInd) = mean(temp);
        sigQErr(indd,snrInd) = std(temp);

        temp = errUloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHU(:,indd,snrInd) = hh/sum(hh);
        meanUErr(indd,snrInd) = mean(temp);
        sigUErr(indd,snrInd) = std(temp);

        temp = errVloc(inds(indStart:indEnd));
        hh = histc(temp,errSbins);
        HHV(:,indd,snrInd) = hh/sum(hh);
        meanVErr(indd,snrInd) = mean(temp);
        sigVErr(indd,snrInd) = std(temp);

        indStart = indEnd+1;
    end
    snrInd
end

%%
HH = (HHQ + HHU + HHV)/3;
meanErr = (meanQErr + meanUErr + meanVErr)/3;
sigErr = (sigQErr + sigUErr + sigVErr)/3;

figure(2)
clf
snrInd = 10;

subplot(2,3,1)
imagesc(Sinabins(1:end-1),errSbins,HH(:,:,snrInd))
xlabel('Sina')
ylabel('errRet')
title(sprintf('errQUV, SNR = %.1f',SNR(snrInd)))
colorbar


subplot(2,3,2)
imagesc(SNR,Sinabins(1:end-1),meanErr,[-.1,.1])
xlabel('SNR')
ylabel('Sina')
colorbar
title(sprintf('meanErrQUV'))


subplot(2,3,3)
imagesc(SNR,Sinabins(1:end-1),sigErr)
xlabel('SNR')
ylabel('Sina')
colorbar
title(sprintf('stdErrQUV'))


dInd = 40;
subplot(2,3,4)
plot(errSbins,HHQ(:,dInd,snrInd))
hold on
plot(errSbins,HHU(:,dInd,snrInd))
plot(errSbins,HHV(:,dInd,snrInd))
xlabel('Error')
title(sprintf('At SNR = %.2f and D = %.2f-%.2f',SNR(snrInd),Sinabins(dInd),Sinabins(dInd+1)))
legend('Q','U','V')


subplot(2,3,5)
plot(Sinabins(1:end-1),meanErr')
xlabel('Sina')
ylabel('Ret')
title(sprintf('Mean error'))

subplot(2,3,6)
plot(Sinabins(1:end-1),sigErr')
xlabel('Sina')
ylabel('Ret')
title(sprintf('Std error'))
sinaModeling