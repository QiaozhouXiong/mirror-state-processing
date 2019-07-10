function f = filterGen(fwx)
%   Generate a gaussian filter with kernel width fwx

nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
if fwx == 1
    nx = 0;
end
h = exp(-nx.^2);
f = h(:)/sum(h(:));