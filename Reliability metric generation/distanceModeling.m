%%
clear A B
xD = (Dbins(1:end-1)+Dbins(1+1:end))/2;ySNR = SNR;
zsigErr = sigErr;
figure(3);clf;subplot(2,3,1);plot(xD,zsigErr');xlabel('D'); ylabel('Ret');title('Std error');hold on; axis([0 2 0 10])
ft = fittype('1/(a*x+b)');
stdTh = linspace(1,1.1,length(ySNR));
for index = 1:numel(ySNR)
    z = zsigErr(:,index);
    ind = z<stdTh(index);
    StartPoint = [1,2];
    [f, gofl] = fit(xD(ind).',z(ind),ft,'StartPoint',StartPoint);
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
xlabel('D')
ylabel('Ret')
title(sprintf('Std error'));
hold on;
%%
Apval = fitA3; Bpval = fitB3;
xD0 = linspace(0.0,2,50);
SNR0 = linspace(6,35,41);
for index = 1:length(SNR0)
    snr = SNR0(index);
    retstdfit = 1./(polyval(Apval,snr)*xD0+polyval(Bpval,snr));
    plot(xD0,retstdfit,'--')
end
axis([0,2,0,100])
save('RelMetdistance','Apval','Bpval')