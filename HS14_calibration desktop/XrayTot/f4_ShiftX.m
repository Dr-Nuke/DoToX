function [imOut] = f4_ShiftX(im,shift)
% shifts an image in x direction, subpixel

if shift==0 % who calls this function???
    imOut=im;
else


%padd left & right
pad=ceil(abs(shift));
x1=im(1,:);
xe=im(end,:);
im2=vertcat(repmat(x1,pad,1),im,repmat(xe,pad,1));

n=size(im2);
[x,y]=ndgrid(1:n(1),1:n(2));
imOut=interpn(x,y,im2,x-shift,y,'cubic');

%crop back
imOut([1:pad,end-(pad-1):end],:)=[];
end

end

