 function [im_out] = f_ffcorrect(im_in,xdose,ydose)
%% dosed correction

% im_OB = open beam image
% im_in = the image (stack) to be corrected
% im_out = the corrected image (stack)

% in principal this function makes sure that the average value within the
% dose correction window is 1.


if length(size(im_in))==3 %its the image stack...
    dosefactor=(length(xdose)*length(ydose))/sum(sum(im_in(xdose,ydose,:)));
    dosefactorStack=repmat(reshape(dosefactor(:),1,1,length(dosefactor)),size(im_in,1),size(im_in,2));

    im_out=im_in.*dosefactorStack;
elseif length(size(im_in))==2 % only a single image
    dosefactor=(length(xdose)*length(ydose))/sum(sum(im_in(xdose,ydose)));
    im_out=im_in*dosefactor;
end
%end

return


; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

; %10

