function [im_out] = f_spotremove1(im_in,kernel,thresh,mode)
%removes spots

filt=f_MediCustom(im_in,kernel,mode);

diff=im_in-filt;

im_in(diff>threshold*im_in)=filt(diff>threshold*im_in);
im_out=im_in;
end

