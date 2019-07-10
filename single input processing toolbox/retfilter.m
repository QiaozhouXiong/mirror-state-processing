function retmean = retfilter(imIn,width,weight,maskIn)
% compute the compund retardation in the ROI defined by 300um diameter
% circle.
%
% input: imIn, the input data waiting for averaging in two dimension
%        width, the width defines the diameter of the circle
%        weight, the reliability weight for the avearging
%        maskIn, the mask for the avearging, equal to weight by default
%
% output: retmean: the averaged retardatio in a array
% latest version: 10 July

if nargin < 4
    maskIn = weight;
end
maskUQ = (maskIn==0);    
%% 20190602: updates the mask for excluded points where half of the points are not qualified in ROI
imIn = imIn.*weight;
h = fspecial('disk',round(width/2));
h = h(1:width,1:width);
lInds = find(h>0);
subMask = zeros(size(weight));
subMask(1:round(width/2):end-width,1:round(width/2):end-width) = 1;
Mask = logical(maskIn.*subMask);
mn = imIn;

for depthInd = 1:round(width/2):(size(imIn,1)-width)
    
    for latInd = 1:round(width/2):(size(imIn,2)-width)
        
        loc = weight((depthInd-1)+(1:width),(latInd-1)+(1:width));
        Wsum = sum(loc(lInds)); 
        zeroIdx = maskUQ((depthInd-1)+(1:width),(latInd-1)+(1:width));
        if sum(zeroIdx(lInds))>numel(lInds)/2
            Wsum = 0;
        end
        loc = imIn((depthInd-1)+(1:width),(latInd-1)+(1:width));
        mn(depthInd,latInd) = sum(loc(lInds))/Wsum;
    end
end
Mask = and(and(Mask,~isnan(mn)),~isinf(mn));
retmean = mn(Mask);
% out.med = med(~isnan(med));
% out.mean = mn(~isnan(med));
% out.max = mx(~isnan(med));