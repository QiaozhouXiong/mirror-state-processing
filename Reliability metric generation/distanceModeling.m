%%
clear A B
D_select_idx = 25;
SNR_select_idx = 10;
xD = (Dbins(D_select_idx:end-1)+Dbins(D_select_idx+1:end))/2;ySNR = SNR(SNR_select_idx:end);
zsigErr = sigErr(D_select_idx:end,SNR_select_idx:end);
figure(11);clf;subplot(2,3,1);plot(xD,zsigErr');xlabel('D'); ylabel('Ret');title(sprintf('Std error'));hold on;
for index = 1:numel(ySNR)
    ft = fittype('1/(a*x+b)');
    StartPoint = [1,2];
    z = zsigErr(:,index);
    [f, gofl] = fit(xD.',z,ft,'StartPoint',StartPoint);
    A(index) = f.a; B(index) = f.b;
%     C(index) = f.c; D(index) = f.d;
    plot(f,'--');legend off
end
subplot(2,3,2);
plot(ySNR,A.');hold on;
xlabel('SNR')
ylabel('fitted a');
fitA3 = polyfit(ySNR,A,3);
plot(ySNR,polyval(fitA3,ySNR))

% B modeling
subplot(2,3,3);
plot(ySNR,B);hold on
plot(ySNR,repmat(mean(B),size(ySNR)))
fitB2 = polyfit(ySNR,B,2);
plot(ySNR,polyval(fitB2,ySNR));
fitB3 = polyfit(ySNR,B,3);
plot(ySNR,polyval(fitB3,ySNR));
legend('b','mean value','2-order fitting','3-order fitting')

%%
subplot(2,3,[4 5 6]);
D_select_idx = 1;
xD_orig = (Dbins(D_select_idx:end-1)+Dbins(D_select_idx+1:end))/2;
plot(xD_orig,sigErr(D_select_idx:end,:))
xlabel('D')
ylabel('Ret')
title(sprintf('Std error'));
hold on;
%%
Apval = fitA3; Bpval = fitB3;
for index = 1:length(SNR)
    snr = SNR(index);
    retstdfit = 1./(polyval(Apval,snr)*xD_orig+polyval(Bpval,snr));
    plot(xD_orig,retstdfit,'--')
end
save('RelMetdistance','Apval','Bpval')