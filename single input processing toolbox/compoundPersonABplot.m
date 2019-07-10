function out = compoundPersonABplot(Ret1, Ret2, Ret3, Rel2, Rel3, relTh, n)
% COmpare the single input processing with two input processing, and
% compare the single input processing in two channels in the way of Perason
% correlation and Bland-Altman plot
% 
% features added in this version
% 1. using scatting to add hue to the measuremnt in regresssion
% 2. add the comparison between two channel single input processing
%
% input args:
% Ret1, Ret2, Ret3, two input procesing, single input processing with s1
% and S2, respectively.
%
% Rel2, Rel3, the averaged measurement reliability of the measurement
%
% relTh, relTh threshold, 0 percent by default
%
% n, figure number, if 0 means generate seperate figures for latter
% processing. single figure by default, namely n = 4.
%
% lastest version updated at 10 July 10 2019 by Qiaozhou Xiong

if nargin<4
    Rel2 = ones(size(Ret2));
    Rel3 = ones(size(Ret3));
end
if nargin <6
    relTh = 0.0;
end
if nargin<7
    n = 4;
end
flag = false;
if n == 0
    flag = true;
    close all;
    n = 1;
end
%% set reliability cut off

a = .95;
figure(n); clf
% calculate the cutoff limit
% Rel2 = mat2gray((1./(sqrt(Rel2))), [0.02, 0.2]);
% Rel3 = mat2gray((1./(sqrt(Rel3))), [0.02, 0.2]);
Rel2 = (1./(sqrt(Rel2)));
Rel3 = (1./(sqrt(Rel3)));
sampleNo = length(Rel3);
Rel2Sorted = sort(Rel2);
Rel3Sorted = sort(Rel3);

Rel2Th = Rel2Sorted(round((1-relTh)*sampleNo));
Rel3Th = Rel3Sorted(round((1-relTh)*sampleNo));

Rel2mask = Rel2<=Rel2Th;
Rel3mask = Rel3<=Rel3Th;
%% compare how similar
if flag
    ax1 = figure(1);
    set(ax1, 'position', [800 400 300 300]);
else
    ax1 = subplot(3,2,1);
end

hold on;
scatter(Ret1,Ret2,2,mat2gray(Rel2,double([prctile(Rel2,10), prctile(Rel2,90)])))
colormap(ax1,ametrine(256)); 
hc = colorbar; 
set(hc,'Ytick', [0 1])
b = deming(Ret1(Rel2mask),Ret2(Rel2mask));
x = linspace(0,2e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(Ret1(Rel2mask),Ret2(Rel2mask));
correlation = corrIdx(2,1);
xlabel('\Deltan_1_2'); ylabel('\Deltan_1'); axis([0,2e-3,0,2e-3])
text(0.2e-3,1.9e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,1.7e-3,sprintf('r=%1.2g',correlation),'Color','black')
%title('ret1 ret2 correlation'); 
axis xy; axis square; box on

if flag
    ax2 = figure(2);
    set(ax2, 'position', [400 400 500 250]);
else
    ax2 = subplot(3,2,2);
end
retbinsx = linspace(0,2e-3,201);
retbinsy = linspace(-.5e-3,.5e-3,101);
hh2 = hist2D((Ret1(Rel2mask)+Ret2(Rel2mask))./2,Ret1(Rel2mask)-Ret2(Rel2mask),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy; colormap(ax2, parula);%axis square
xlabel('(\Deltan_1_2 +\Deltan_1)/2');ylabel('\Deltan_1_2-\Deltan_1');
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((Ret1(Rel2mask)-Ret2(Rel2mask)));
retmean = sort((Ret1(Rel2mask)+Ret2(Rel2mask))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
plot(retAvg2*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
text(double(retmeanMed-1e-4),-0.45e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;
plot(linspace(0,2e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',1);
text(1.2e-3,double(retmeanDiff+1e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2e-3,500),LoAplus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAplus+1e-4),'LoA(+)','Color','w')
text(1.2e-3,double(LoAplus+1e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2e-3,500),LoAminus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAminus-.5e-4),'LoA(-)','Color','w')
text(1.2e-3,double(LoAminus-.5e-4),num2str(LoAminus,2),'Color','w');
%%
out.slope1 = b(2);
out.Corr1 = correlation;
out.LoA1 = LoAplus-LoAminus;
%%

%% compare how similar
if flag
    ax3 = figure(3);
    set(ax3, 'position', [800 400 300 300]);
else
    ax3 = subplot(3,2,3);
end
hold on;
scatter(Ret1,Ret3,2,mat2gray(Rel3,double([prctile(Rel3,10), prctile(Rel3,90)])))
colormap(ax3, ametrine(256));
hc = colorbar; set(hc,'Ytick', [0 1])
b = deming(Ret1(Rel3mask),Ret3(Rel3mask));
x = linspace(0,2e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(Ret1(Rel3mask),Ret3(Rel3mask));
correlation = corrIdx(2,1);
xlabel('\Deltan_1_2'); ylabel('\Deltan_2'); axis([0,2e-3,0,2e-3])
text(0.2e-3,1.9e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,1.7e-3,sprintf('r=%1.2g',correlation),'Color','black')
% title('ret1 ret3 correlation');
axis xy; axis square; box on

if flag
    ax4 = figure(4);
    set(ax4, 'position', [400 400 500 250]);
else
    ax4 = subplot(3,2,4);
end
retbinsx = linspace(0,2e-3,201);
retbinsy = linspace(-.5e-3,.5e-3,101);
hh2 = hist2D((Ret1(Rel3mask)+Ret3(Rel3mask))./2,Ret1(Rel3mask)-Ret3(Rel3mask),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy;colormap(ax4,parula) %axis square
xlabel('(\Deltan_1_2 +\Deltan_2)/2');ylabel('\Deltan_1_2-\Deltan_2');
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((Ret1(Rel3mask)-Ret3(Rel3mask)));
retmean = sort((Ret1(Rel3mask)+Ret3(Rel3mask))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
plot(retAvg2*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
text(double(retmeanMed-1e-4),-0.45e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;

plot(linspace(0,2e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',1);
text(1.2e-3,double(retmeanDiff+1e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2e-3,500),LoAplus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAplus+1e-4),'LoA(+)','Color','w')
text(1.2e-3,double(LoAplus+1e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2e-3,500),LoAminus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAminus-.5e-4),'LoA(-)','Color','w')
text(1.2e-3,double(LoAminus-.5e-4),num2str(LoAminus,2),'Color','w');
out.slope2 = b(2);
out.Corr2 = correlation;
out.LoA2 = LoAplus-LoAminus;
%% compare how similar
if flag
    ax5 = figure(5);
    set(ax5, 'position', [800 400 300 300]);
else
    ax5 = subplot(3,2,5);
end

hold on;
scatter(Ret2,Ret3,2,mat2gray((Rel2+Rel3)./2,double([prctile((Rel2+Rel3)./2,10), prctile((Rel2+Rel3)./2,90)])))
colormap(ax5, ametrine(256));
hc = colorbar; set(hc,'Ytick', [0 1])
b = deming(Ret2(Rel3mask),Ret3(Rel3mask));
x = linspace(0,2e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(Ret2(Rel3mask),Ret3(Rel3mask));
correlation = corrIdx(2,1);
xlabel('\Deltan_1'); ylabel('\Deltan_2'); axis([0,2e-3,0,2e-3])
text(0.2e-3,1.9e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,1.7e-3,sprintf('r=%1.2g',correlation),'Color','black')
% title('ret1 ret3 correlation');
axis xy; axis square; box on

if flag
    ax6 = figure(6);
    set(ax6, 'position', [400 400 500 250]);
else
    ax6 = subplot(3,2,6);
end
retbinsx = linspace(0,2e-3,201);
retbinsy = linspace(-.5e-3,.5e-3,101);
hh2 = hist2D((Ret2(Rel3mask)+Ret3(Rel3mask))./2,Ret2(Rel3mask)-Ret3(Rel3mask),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy;colormap(ax4,parula) %axis square
xlabel('(\Deltan_1 +\Deltan_2)/2');ylabel('\Deltan_1-\Deltan_2');
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((Ret2(Rel3mask)-Ret3(Rel3mask)));
retmean = sort((Ret2(Rel3mask)+Ret3(Rel3mask))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
plot(retAvg2*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
text(double(retmeanMed-1e-4),-0.45e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;

plot(linspace(0,2e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',1);
text(1.2e-3,double(retmeanDiff+1e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2e-3,500),LoAplus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAplus+1e-4),'LoA(+)','Color','w')
text(1.2e-3,double(LoAplus+1e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2e-3,500),LoAminus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAminus-.5e-4),'LoA(-)','Color','w')
text(1.2e-3,double(LoAminus-.5e-4),num2str(LoAminus,2),'Color','w');
out.slope3 = b(2);
out.Corr3 = correlation;
out.LoA3 = LoAplus-LoAminus;