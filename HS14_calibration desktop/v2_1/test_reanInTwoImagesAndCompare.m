
im1=im2double(fitsread('C:\data\tomo_HS14\processed\1raw\ob_\ob_1.fits')');
im1=im1(81:end,201:2450);
im1=im1(6:end-5,6:end-5);
size(im1)


%im1=im1(300:350,600:700);

im2=              load('C:\data\tomo_HS14\processed\2flt\ob_\ob_filt1.mat');
im2=im2.var;
size(im2)
%im2=im2(300:350,600:700);

cmin=0;
cmax=1;
imbo3(im1,1);
caxis([cmin cmax])
imbo3(im2,2);
imbo3((im2-im1)~=0,3);
a=0.1
imbo3((im1*a+im2*(1-a)),4);




; %1
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


