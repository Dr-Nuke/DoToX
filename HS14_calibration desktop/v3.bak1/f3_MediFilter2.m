function [ im_out] = f3_medifilter2( im_in,kernel, threshold,mode )
% median filter for intermediate sized bright and dark spots (up to 8 pixel
% wide)


%apply filter
filt=f_MediCustom(im_in,kernel,mode);

% make the difference
diff=im_in-filt;

% correct output
im_in(abs(diff)>threshold*im_in)=filt(abs(diff)>threshold*im_in);

im_out=im_in;

end