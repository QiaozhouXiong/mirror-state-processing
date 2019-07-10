function out = comVisualization(int,out2input,out1input,pst,n)
% Visualize the comparison between the benchmark two input processing(out2input) and single input processing(out1input).
%
% input:  int: intensity image of the slice
%         pst: ps processing parameters
%         pst.offset: the offset to tranform the figure into cartisian
%         coordinate
%         n: the figure index, 2 by default. if n is eual to n, generate
%         the sepearte figures for latter processing
%
% output: The compound retardation in the ROI defined by 300um diameter of dual inputs and single input,
%         and compound reliability of single input.
% 
%%-------------------------------------------------------------------------

if nargin<5
    n = 2;
end
offset = pst.offset;
clipLimit = pst.clipLimit;
figure(n)
clf;
%% read the data for comparison
ret1 = out2input.rmeancorr2; % conventional two input processing
ret2 = out1input.ret1;        % single input processing
ret3 = out1input.ret2;        % single input processing
dop = out2input.dop;
%% define  the area for comparison
mask = dop>.8;
mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;
%% 
logLim = [55,110];
im1 = signal2isolum(ret1,10*log10(int),[0,0.0024],logLim);
im2 = signal2isolum(ret2,10*log10(int),[0,0.0024],logLim);
im3 = signal2isolum(ret3,10*log10(int),[0,0.0024],logLim);
Nout = 512;
Range = offset+ [0, 600];
cim1 = circInterp(im1,Nout,Range);
cim2 = circInterp(im2,Nout,Range);
cim3 = circInterp(im3,Nout,Range);
cdop = circInterp(dop,Nout,Range);
cint = min(max((circInterp(10*log10(int),Nout,Range)-logLim(1))/diff(logLim),0),1);
out.cint = cint;

nonlinMapping = [170,20];
temp = exp(-((0:255)-nonlinMapping(1)).^2/nonlinMapping(2)^2);
temp = cumsum(temp)/sum(temp);% erf
mask1 = repmat(temp(round(cdop*255)+1),[1,1,3]);
mask2 = 1-mask1;

out.cim1m = cim1.*mask1 + repmat(cint,[1,1,3]).*mask2;
out.cim2m = cim2.*mask1 + repmat(cint,[1,1,3]).*mask2;
out.cim3m = cim3.*mask1 + repmat(cint,[1,1,3]).*mask2;
% to determin the area to compare
ret1 = ret1.*mask;
ret2 = ret2.*mask;
ret3 = ret3.*mask;
ret3(isnan(ret3)) = 0;
ret2(isnan(ret2)) = 0;
ret1(isnan(ret1)) = 0;
%% reliability metric
%convert reliablity to estimated std
Rel = sqrt((sum(out1input.avgWeg1,3)));
rel2 = mat2gray(sqrt(1./Rel),double([mean(prctile(sqrt(1./Rel),10)) mean(prctile(sqrt(1./Rel),90))]));
Rel = sqrt(sum(out1input.avgWeg2,3));
rel3 = mat2gray(sqrt(1./Rel),double([mean(prctile(sqrt(1./Rel),10)) mean(prctile(sqrt(1./Rel),90))]));
rel23 = (rel2+rel3)/2;
%% can set a cut-off reliability threshold here
mask2 = and(mask,(rel2>=0.0));
mask3 = and(mask,(rel3>=0.0));
mask23 = and(mask2,mask3);
%% compare how similar/channle 1 and S12
if n==1
    ax1 = subplot(2,3,1);
else
    ax1 = figure(1);set(gcf,'position',[200 200 400 200])
end
hold on;
scatter(ret1(mask2),ret2(mask2),1,rel2(mask2))
colormap(ax1,ametrine(256)); colorbar;caxis([0 1])
b = deming(ret1(mask2),ret2(mask2));
x = linspace(0,2.5e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
axis([0 2.5e-3, 0 2.5e-3])
corrIdx = corrcoef(ret1(mask2),ret2(mask2));
correlation = corrIdx(2,1);
xlabel('\Deltan_1_2'); ylabel('\Deltan_1');
text(0.2e-3,2.4e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,2.2e-3,sprintf('r=%1.2g',correlation),'Color','black')
axis xy; axis square
out.correlation(1) = correlation;
out.slope(1) = b(2);
%% channle 2 and S12
if n==1
    ax2 = subplot(2,3,2);
else
    ax2 = figure(2);set(gcf,'position',[200 200 400 200])
end
hold on;
scatter(ret1(mask3),ret3(mask3),1,rel3(mask3))
colormap(ax2,ametrine(256)); colorbar; caxis([0 1])
b = deming(ret1(mask3),ret3(mask3));
x = linspace(0,2.5e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(ret1(mask3),ret3(mask3));
correlation = corrIdx(2,1);
xlabel('\Deltan_1_2'); ylabel('\Deltan_2'); axis([0,2.5e-3,0,2.5e-3])
text(0.2e-3,2.4e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,2.2e-3,sprintf('r=%1.2g',correlation),'Color','black')
axis xy; axis square
out.correlation(2) = correlation;
out.slope(2) = b(2);
%% channle 1 and channel 2
if n==1
    ax3 = subplot(2,3,3);
else
    ax3 = figure(3);set(gcf,'position',[200 200 400 200])
end

hold on;
scatter(ret2(mask23),ret3(mask23),1,rel23(mask23))
colormap(ax3,ametrine(256)); colorbar; caxis([0 1])
b = deming(ret2(mask23),ret3(mask23));
x = linspace(0,2.5e-3,200);
y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(ret2(mask23),ret3(mask23));
correlation = corrIdx(2,1);
xlabel('\Deltan_1'); ylabel('\Deltan_2'); axis([0,2.5e-3,0,2.5e-3])
text(0.2e-3,2.4e-3,sprintf('a=%1.2g',b(2)),'Color','black')
text(0.2e-3,2.2e-3,sprintf('r=%1.2g',correlation),'Color','black')
axis xy; axis square
out.correlation(3) = correlation;
out.slope(3) = b(2);
%% Bland-Atman Analysis
a = .95;
if n==1
    ax4 = subplot(2,3,4);
else
    ax4 = figure(4);set(gcf,'position',[200 200 400 200])
end
retbinsx = linspace(0,2.5e-3,201);
retbinsy = linspace(-1.5e-3,1.5e-3,201);
hh2 = hist2D((ret1(mask2)+ret2(mask2))./2,ret1(mask2)-ret2(mask2),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy; colormap(ax4, parula);%axis square
xlabel('(\Deltan_1_2 +\Deltan_1)/2');ylabel('\Deltan_1_2-\Deltan_1');
set(gca,'Xtick',[0 0.5 1 1.5 2].*1e-3);set(gca,'Ytick',[-1 -0.5 0 0.5 1].*1e-3);axis([0 2.5 -1.5 1.5]*1e-3)
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((ret1(mask2)-ret2(mask2)));
retmean = sort((ret1(mask2)+ret2(mask2))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
plot(retAvg2*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
text(double(retmeanMed+1e-4),-0.65e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;
plot(linspace(0,2.5e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2.5e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',.5);
text(2e-3,double(retmeanDiff+2e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAplus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAplus+2e-4),'LoA(+)','Color','w')
text(2e-3,double(LoAplus+2e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAminus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAminus-1.5e-4),'LoA(-)','Color','w')
text(2e-3,double(LoAminus-1.5e-4),num2str(LoAminus,2),'Color','w');
%%
if n==1
    ax5 = subplot(2,3,5);
else
    ax5 = figure(5);set(gcf,'position',[200 200 400 200])
end
retbinsx = linspace(0,2.5e-3,201);
retbinsy = linspace(-1.5e-3,1.5e-3,201);
hh2 = hist2D((ret1(mask3)+ret3(mask3))./2,ret1(mask3)-ret3(mask3),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy; colormap(ax5, parula);%axis square
xlabel('(\Deltan_1_2 +\Deltan_2)/2');ylabel('\Deltan_1_2-\Deltan_2');
set(gca,'Xtick',[0 0.5 1 1.5 2].*1e-3);set(gca,'Ytick',[-1 -0.5 0 0.5 1].*1e-3);axis([0 2.5 -1.5 1.5]*1e-3)
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((ret1(mask3)-ret3(mask3)));
retmean = sort((ret1(mask3)+ret3(mask3))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
plot(retAvg2*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
text(double(retmeanMed+1e-4),-0.65e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;
plot(linspace(0,2.5e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2.5e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',.5);
text(2e-3,double(retmeanDiff+2e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAplus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAplus+2e-4),'LoA(+)','Color','w')
text(2e-3,double(LoAplus+2e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAminus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAminus-1.5e-4),'LoA(-)','Color','w')
text(2e-3,double(LoAminus-1.5e-4),num2str(LoAminus,2),'Color','w');
%%
if n==1
    ax6 = subplot(2,3,6);
else
    ax6 = figure(6);
    set(gcf,'position',[200 200 400 200])
end
retbinsx = linspace(0,2.5e-3,201);
retbinsy = linspace(-1.5e-3,1.5e-3,201);
hh2 = hist2D((ret2(mask23)+ret3(mask23))./2,ret2(mask23)-ret3(mask23),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log10(hh2')); hold on;axis xy; colormap(ax6, parula);%axis square
xlabel('(\Deltan_1 +\Deltan_2)/2');ylabel('\Deltan_1-\Deltan_2');
set(gca,'Xtick',[0 0.5 1 1.5 2].*1e-3);set(gca,'Ytick',[-1 -0.5 0 0.5 1].*1e-3);axis([0 2.5 -1.5 1.5]*1e-3)
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((ret2(mask23)-ret3(mask23)));
retmean = sort((ret2(mask23)+ret3(mask23))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
plot(retAvg2*ones(500,1),linspace(-1.5e-3,1.5e-3,500),'w--','lineWidth',.5);
text(double(retmeanMed+1e-4),-0.85e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;
plot(linspace(0,2.5e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2.5e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',.5);
text(2e-3,double(retmeanDiff+2e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAplus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAplus+2e-4),'LoA(+)','Color','w')
text(2e-3,double(LoAplus+2e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAminus*ones(500,1),'w--','lineWidth',.5);
text(0.0,double(LoAminus-1.5e-4),'LoA(-)','Color','w')
text(2e-3,double(LoAminus-1.5e-4),num2str(LoAminus,2),'Color','w');
if n==1
    set(gcf,'position',[200 200 1500 400])
end 
%% calculation
ret1 = circInterp(ret1,Nout,Range);
ret2 = circInterp(ret2,Nout,Range);
ret3 = circInterp(ret3,Nout,Range);
maskC = circInterp(mask,Nout,Range);
width = 34;
rel2C = circInterp(sum(out1input.avgWeg1,3),Nout,Range);
rel3C = circInterp(sum(out1input.avgWeg2,3),Nout,Range);

rel2Cs = retfilter(rel2C,width,maskC);
rel3Cs = retfilter(rel3C,width,maskC);

% rel2mask = mat2gray((rel2),[0,10]);
% rel3mask = mat2gray((rel3),[0,10]);
% h = fspecial('disk',width/2);
% h = h(1:width,1:width);
ret1s = retfilter(ret1,width,maskC);
ret2s = retfilter(ret2,width,(maskC.*rel2C),maskC);
ret3s = retfilter(ret3,width,(maskC.*rel3C),maskC);

out.rel2 = rel2Cs;
out.rel3 = rel3Cs;
out.ret1 = ret1s;
out.ret2 = ret2s;
out.ret3 = ret3s;