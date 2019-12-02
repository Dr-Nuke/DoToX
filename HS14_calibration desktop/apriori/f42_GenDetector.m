function [ output_args ] = f42_GenDetector(xy0,phi,res,npix)
% generates the detector pixel coordinates

n=10; %number of pixels
res=1; % pixel spacing in mm
l=n*res;       % pixel row length


x=zeros(1,n);
y=linspace(-(n-1)/2*res,(n-1)/2*res,n);
xy=[x;y]';


end

