% filter development
close all

set_folder='chcl3\';
set_file='CHCL3_';
readpath='C:\data\tomo_HS14\02_rawdata\';
writepath='C:\data\tomo_HS14\processed\';


i=1500;
fnn=sprintf('%04d',i); %file name numberfi
name_fits=strcat(readpath,set_folder,set_file,fnn,'.fits');

im=im2double(fitsread(name_fits)');

im2=f_CamLineFilt(im);

im3=f_medifilter1(im2,thresh_medi);

im4=f_MediCustom(im3,thresh_medi,'mean');




x1=600;x2=1200;y1=300,y2=1000;
figure('Name','small spots','units','normalized','outerposition',[0 0 1 1])
imbo4(im3');axis([x1,x2,y1,y2]);axis equal;


figure('Name','Camlines','units','normalized','outerposition',[0 0 1 1])
imbo4(im2');axis([x1,x2,y1,y2]);;axis equal;


fig=figure('Name','Raw image','units','normalized','outerposition',[0 0 1 1])
imbo4(im');axis([x1,x2,y1,y2]);;axis equal;

% fig=figure('Name','DRaw','units','normalized','outerposition',[0 0 1 1])
% imbo4((im2.*(im3~=im2))');axis([x1,x2,y1,y2]);;axis equal;


figure('Name','customSpot','units','normalized','outerposition',[0 0 1 1])
imbo4(im4');axis([x1,x2,y1,y2]);axis equal;