 function [im_out] = f_ffcorrect(im_OB,im_in)
%% flat field correction
% im_OB = open beam image
% im_in = the image (stack) to be corrected
% im_out = the corrected image (stack

if length(size(im_in))==3 %its the image stack...

    tmp     = repmat(im_OB,1,1,size(im_in,3)); % correction stack
    im_out  = im_in./tmp; % devide image stack by correction stack
    
elseif length(size(im_in))==2 % ...or only an image
    im_out=im_in./im_OB;
end
return


