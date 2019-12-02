
% finds the number of angles
%% Open Files
clear all
close all
clc


filepath=('C:/data/20161216 Campaign 5 Robert Z/11h34 Tmax 110kV 5mA.seq');
fileInfo = dir(filepath);
imsize=[1024 640];
header_size=2048;
fileSize = fileInfo.bytes;
fileSize=floor((fileSize-header_size)/imsize(1)/imsize(2)/2);
fid=fopen(filepath);
hed = fread(fid,header_size);%the header size migth need to be adjusted depending on image settings
%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
img=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16')));
fclose(fid);

%rotate to real geometry
img=rot90(reshape(img,[1024 640 fileSize]));


