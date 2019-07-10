function out = twoinputCom(out2input,out2inputD,pst,n)
% compare the results of subsampled two input processing with benchmark
% processing.
%
% inputs: out2input, the results of dual input processing
%         out2inputD, the subsampled couterparts
%         pst. processing parameters
%         n, figure index, n = 1 by default.
%
% output: the compound averaged retardation
%
% latest version at 10 July
%--------------------------------------------------------------------------

if nargin<4
    n = 1;
end
offset = pst.offset;
clipLimit = pst.clipLimit;
figure(n)
clf
ret0 = out2input.rmeancorr2; % conventional two input processing
ret1 = out2inputD.rmeancorr2; % conventional two input processing

mask = out2input.dop>.8;
mask(1:clipLimit,:) = 0;
mask(end-clipLimit:end,:) = 0;
% h = [0.5 0.5];
% retf = imfilter(ret0,h,'replicate');
ret0 = ret0(:,1:2:end);
mask = mask(:,1:2:end);
% offset = 0;
Range = offset+ [0, 800];
% to determin the area to compare
ret0 = ret0.*mask;
ret1 = ret1.*mask;
ret0(isnan(ret0)) = 0;
ret1(isnan(ret1)) = 0;
%% deming correlation
retbinsx = linspace(0,2e-3,201);
retbinsy = linspace(0,2e-3,201);
ax1 = subplot(1,2,1); hold on;
hh1 = hist2D(ret0(mask),ret1(mask),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log(hh1'));colormap(ax1,ametrine(256));
b = deming(ret0(mask),ret1(mask));
x = linspace(0,2e-3,200);y = x*b(2)+b(1);
plot(x,y,'k-','linewidth',0.5);
corrIdx = corrcoef(ret0(mask),ret1(mask));
correlation = corrIdx(2,1);
xlabel('\Deltan1'); ylabel('\Deltan2'); axis([0,2e-3,0,2e-3])
text(0.2e-3,1.9e-3,sprintf('a=%1.2g',b(2)),'Color','green')
text(0.2e-3,1.7e-3,sprintf('r=%1.2g',correlation),'Color','green')
axis xy; axis square
% title('\Deltan1 2 correlation'); 
out.correlation(1) = correlation;
out.slope = b(2);
%% Bland-Atman Analysis
a = .95;
%%
retbinsx = linspace(0,2.5e-3,201);
retbinsy = linspace(-0.8e-3,0.8e-3,201);
ax2 = subplot(1,2,2);
hh2 = hist2D((ret0(mask)+ret1(mask))./2,ret0(mask)-ret1(mask),retbinsx,retbinsy);
imagesc(retbinsx,retbinsy,log(hh2')); hold on;axis xy; axis equal; colormap(ax2, parula);%axis square
xlabel('(\Deltan_1 +\Deltan_2)/2');ylabel('\Deltan_1-\Deltan_2');
set(gca,'Xtick',[0 0.5 1 1.5 2].*1e-3);set(gca,'Ytick',[-1 -0.5 0 0.5 1].*1e-3);axis([0 2.5 -0.8 0.8]*1e-3)
hc = colorbar; set(hc,'Ytick', [])
% calculation
retDiff = sort((ret0(mask)-ret1(mask)));
retmean = sort((ret0(mask)+ret1(mask))./2);
sampleNo = length(retmean);

retAvg1 = retmean(round((1-a)/2*sampleNo));
retAvg2 = retmean(round((1+a)/2*sampleNo));
retmeanMed = retmean(round(1/2*sampleNo));
retmeanAvg = median(retmean)
plot(retAvg1*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
plot(retAvg2*ones(500,1),linspace(-1e-3,1e-3,500),'w--','lineWidth',1);
text(double(retmeanMed+1e-4),-0.65e-3,num2str(retAvg2-retAvg1,2),'Color','w')

retmeanDiff = retDiff(round(1/2*sampleNo));
LoAplus = retDiff(round((1+a)/2*sampleNo))-retmeanDiff;
LoAminus = retDiff(round((1-a)/2*sampleNo))-retmeanDiff;
plot(linspace(0,2.5e-3,500),zeros(500,1),'w-','lineWidth',1);
plot(linspace(0,2.5e-3,500),retmeanDiff*ones(500,1),'w--','lineWidth',1);
text(2e-3,double(retmeanDiff+1e-4),num2str(retmeanDiff,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAplus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAplus+1e-4),'LoA(+)','Color','w')
text(2e-3,double(LoAplus+1e-4),num2str(LoAplus,2),'Color','w')
plot(linspace(0,2.5e-3,500),LoAminus*ones(500,1),'w--','lineWidth',1);
text(0.0,double(LoAminus-.5e-4),'LoA(-)','Color','w')
text(2e-3,double(LoAminus-.5e-4),num2str(LoAminus,2),'Color','w');
set(gcf,'position',[200 200 1000 250])

%% calculation
Nout = 512;
ret1 = circInterp(ret1,Nout,Range);
ret0 = circInterp(ret0,Nout,Range);
maskC = circInterp(mask,Nout,Range);
width = 17;
% h = fspecial('disk',width/2);
% h = h(1:width,1:width);
ret1s = retfilter(ret1,width,maskC);
ret0s = retfilter(ret0,width,maskC);
out.ret1cat = ret1;
out.ret0cat = ret0;
out.ret1 = ret1s;
out.ret0 = ret0s;