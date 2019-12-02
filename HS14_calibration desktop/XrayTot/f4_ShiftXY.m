function [imOut] = f4_ShiftXY(im,shift)
% shifts an image in x direction, subpixel
% im = image 
% shift = [xshift,yshift]


%padd left & right
xpad=ceil(abs(shift(1)));
ypad=ceil(abs(shift(2)));

x1=im(1,:);
xe=im(end,:);
im2=vertcat(repmat(x1,xpad,1),im,repmat(xe,xpad,1));

y1=im2(:,1);
ye=im2(:,end);
im3=horzcat(repmat(y1,1,ypad),im2,repmat(ye,1,ypad));


n=size(im3);
[x,y]=ndgrid(1:n(1),1:n(2));
imOut=interpn(x,y,im3,x-double(shift(1)),...
                      y-double(shift(2)),'cubic');

%crop back
  imOut([1:xpad,end-(xpad-1):end],:)=[];
imOut(:,[1:ypad,end-(ypad-1):end])=[];


end
