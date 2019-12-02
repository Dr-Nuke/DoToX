% this script does the x-ray analysis for one file


% this is the 3rd version of total,
% for the slice wise processing and intermediate saving to hard disk
% XRAY Data

clc             % clear commad line output
clearvars       % clear variables
clear global    % clear globals
format compact  % compact command line output

%get name of thisscript
[~,name,~]=fileparts(mfilename('fullpath')); %get m file name

% add path so i can use the same functions
% win
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3_Xray2')
scriptpath='C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3_Xray2/';
dir_orig ='C:/data/20161216 Campaign 5 Robert Z/';
dir_work ='C:/data/20161216 Campaign 5 Robert Z/processed/';

% redhat
% addpath('/home/cbolesch/HS14_calibration desktop')
% addpath('/home/cbolesch/HS14_calibration desktop/v2_1')
% addpath('/home/cbolesch/HS14_calibration desktop/v3')
% scriptpath='/home/cbolesch/HS14_calibration desktop/v3/';
% dir_orig ='/media/data/cbolesch/Tomo_HS14/02_rawdata/';
% dir_work ='/media/data/cbolesch/Tomo_HS14/processed/';

% pre- settings
t_tot=tic;
tic;


%% Open Files

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
img=rot90(reshape(img,[1024 640 fileSize]));




