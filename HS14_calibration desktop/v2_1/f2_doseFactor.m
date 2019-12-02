function [df] = f2_doseFactor(im,xrange)
%S2 Summary of this function goes here
%   Detailed explanation goes here
    y=size(im,2);
    df=sum(sum(im(xrange,:)))/(length(xrange)*y);
    
end

