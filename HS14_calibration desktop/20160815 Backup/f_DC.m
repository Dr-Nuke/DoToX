function [out] = f_DC(im,DC)
%substracts the darc current or scatter contribution
if length(size(im))==3 %its the image stack...
    DCstack=repmat(DC,1,1,size(im,3));
    out=im-DCstack;
    % less or equal zero correction
    
elseif length(size(im))==2
    out=im-DC;
end
out(out<=0)=0.001;

end

