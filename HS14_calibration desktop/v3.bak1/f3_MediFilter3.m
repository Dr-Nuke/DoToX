function [ im_out] = f3_MediFilter3( im_in, threshold, kernel )
% performs the median filter to remove 1-pixel gamma spots
% corrects horizontal lines

%apply filter
filt=medfilt2(im_in,kernel);

% make the difference
diff=im_in-filt;

%copy image
im_out=im_in;

% correct output
im_out(diff>threshold*im_in)=filt(diff>threshold*im_in);



end

