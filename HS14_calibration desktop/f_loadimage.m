function [im] = f_loadimage(i)
%F_LOADIMAGE Summary of this function goes here
set_folder='chcl3\';
set_file='CHCL3_';
readpath='C:\data\tomo_HS14\02_rawdata\';
writepath='C:\data\tomo_HS14\processed\';


i=1500;
fnn=sprintf('%04d',i); %file name numberfi
name_fits=strcat(readpath,set_folder,set_file,fnn,'.fits');

im=im2double(fitsread(name_fits)');
end

