function [ im_out] = f_medifilter1_2( im_in, threshold )
% performs the median filter to remove 1-pixel gamma spots
% corrects vertical lines

%apply filter
filt=medfilt2(im_in,[3 1]);

% make the difference
diff=im_in-filt;

% correct output
im_in(diff>threshold*im_in)=filt(diff>threshold*im_in);

im_out=im_in;

end