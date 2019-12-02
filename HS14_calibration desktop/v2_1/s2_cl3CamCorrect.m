%% 
% his script removes the camera artifact starting at cl3 image 1552
% 
% 
% load('xls_file')
% flag={'emp';'d2o';'cl3'};
% i=3;
% for j=1:

% cl3=fitsread('C:\data\tomo_HS14\02_rawdata\chcl3\CHCl3_0001.fits')';
% cl3_360=fitsread('C:\data\tomo_HS14\02_rawdata\chcl3\CHCl3_1611.fits')';
% ob_cl3=fitsread('C:\data\tomo_HS14\processed\1raw\cl3\cl3_1619.fits')';
% ob_d2o=fitsread('C:\data\tomo_HS14\processed\1raw\d2o\d2o_1619.fits')';
% ob_emp=fitsread('C:\data\tomo_HS14\processed\1raw\emp\emp_0406.fits')';

ob_emp=f2_boLoad('C:\data\tomo_HS14\processed\3add\emp\emp_add_OB__.mat');
ob_d2o=f2_boLoad('C:\data\tomo_HS14\processed\3add\d2o\d2o_add_OB__.mat');
cl3=f2_boLoad('C:\data\tomo_HS14\processed\2flt\cl3\cl3_filt_1583.mat');
ob_cl3=f2_boLoad('C:\data\tomo_HS14\processed\3add\cl3\cl3_add_OB__.mat');


quot=ob_d2o./ob_cl3;
quot2=ob_d2o./cl3;
%%
figure(10);cla; hold on;

m1=mean(ob_cl3./ob_d2o,2);

m2=mean(ob_d2o./ob_emp,2);
m3=mean(ob_cl3./ob_emp,2);



plot(m1,'b');
plot(m2,'r');
plot(m3,'g');

camcorr=ob_d2o./ob_cl3;
cl3_corr=cl3.*camcorr;

%%
caxmin=0.5
imbo3(cl3,1);colorbar
imbo3(cl3_corr,2);colorbar
imbo3(1./quot,10);colorbar;caxis([caxmin 1.2])
imbo3(1./quot2,11);colorbar;caxis([caxmin 1.2])
imbo3(cl3./ob_cl3,12);colorbar;caxis([caxmin 1.2])
;;


%the middle is between pixels 601 and 602
% divide 

; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

;%% %10;;