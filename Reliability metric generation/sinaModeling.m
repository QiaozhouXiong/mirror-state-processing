figure(12);clf
subplot(2,3,1)
plot(Sinabins(1:end-1),sigErr)
xlabel('Sina')
ylabel('Ret')
title(sprintf('Std error'));
hold on; legend off
for index = 1:numel(sigErr(1,:))
    ft = fittype('1/(a*x+b)');
    StartPoint = [1,1];
    x = Sinabins(1:end-1);
    y = sigErr(:,index);
    [f, gofl] = fit(x.',y,ft,'StartPoint',StartPoint);
    A(index) = f.a; B(index) = f.b;
    plot(f,'--'); legend off
end
subplot(2,3,2);
plot(SNR,A);hold on;
xlabel('SNR')
ylabel('fitted a')
x = SNR;
y = A;
% [f,gof] = fit(x.',y.',fittype('a*x+b'),'StartPoint',[1 1]);
p = polyfit(x,y,3);
plot(SNR,polyval(p,SNR))

%%
subplot(2,3,3);
plot(SNR,B);hold on
plot(SNR,repmat(mean(B),size(SNR)))
fitp = polyfit(SNR,B,2);
plot(SNR,polyval(fitp,SNR));
fitp = polyfit(SNR,B,3);
plot(SNR,polyval(fitp,SNR));
legend('b','mean value','2-order fitting','3-order fitting')
% we will choose 3-order fitting here
% B = -1e-4*SNR.^3+3e-3*SNR^2-0.0167*SNR+0.5595;
%% so the equation can be derived as
% f = 1/((1.5608*SNR-0.4193)*sina+0.5383)
% to varify this equation
subplot(2,3,[4 5 6]);
plot(Sinabins(1:end-1),sigErr')
xlabel('Sina')
ylabel('Ret')
title(sprintf('Std error'));
hold on;
x = Sinabins(1:end-1);
for index = 1:length(SNR)
    snr = SNR(index);
    b = fitp(1)*snr.^3+fitp(2)*snr^2+fitp(3)*snr+fitp(4);
%     b = 0.5379
    retstdfit = 1./(polyval(p,snr)*x+polyval(fitp,snr));
    plot(Sinabins(1:end-1),retstdfit,'--')
end
Apval = p;
Bpval = fitp;
save('RelMetsina','Apval','Bpval')
