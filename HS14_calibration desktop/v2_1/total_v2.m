% this is the new version of total,
% for the slice wise processing and intermediate saving to hard disk

clc
clear all
format compact

% add path so i can use the same functions
addpath('V:\1phd\codes\reconstruction\HS14_calibration')


% remarks:
% one image corresponds to 20 seconds of exposure

%% pre- settings

% directories
tic
scriptpath='V:\1phd\codes\reconstruction\HS14_calibration\v2_1\';
readpath ='C:\data\tomo_HS14\02_rawdata\';
dir_orig ='C:\data\tomo_HS14\02_rawdata\';
writepath='C:\data\tomo_HS14\processed\';
dir_work ='C:\data\tomo_HS14\processed\';


%maximum pixel numbers (size of raw files)
xmax=1200; ymax=2450; 
xmin=1; ymin =1;

xpix=[xmin:xmax]; ypix=[ymin:ymax];

%we only need x=80:end since for x<80 the DC image is bad.
%cropping is done in s2_ImageFilter


% threshold for median filter [if (raw-filter) > thresh_medi*raw => correct]
thresh_medi1=0.025; % for gamma spots (single pixel)
thresh_medi2=0.06; % for large spots
kernel_medi2=f_KernelGen(9,9,7.5);
randpix=5; % outer image pixels removed after filtering

% pixel resolution
res=22.1; %pixel/mm, found in resolution.m

% angle list for iradon
angles=(linspace(0,360,376));angles=angles(1:end-1);


% pixel window for dose correction calculation. needs to be outside of any
% structures on all the images

xdo_1 = 300; xdo_2=1050;
xdose=[1:xdo_1-xmin,xdo_2-xmin:xmax-xmin]; 
ydose=[1:length(ypix)];

% write files or just calculate them:
ind_write=1; %

% case separation
itr_cas=[1,2,3,4,5]; % emp, d20, cl3, dc, op
itr_cas=[1];
itr_stp=[1,2,3,4,5]; % ren, flt, add, cor, sin
itr_stp=[1];
dir_raw={'empty\','D2O','chcl3\','CD\','openbeam\'};
pfx_cas={'emp','d2o','cl3','dc_','ob'}; 
dir_cas={'emp\','d2o\','cl3\','dc_\','ob\'};
dir_stp={'1raw\','2flt\','3add\','4cor\','5sin\'};
nfiles={1:406;          % 1:406
        {1:938,1:681};  % {1:938,1:681}  mind special 
        1:1630;         % 1:1630
        1:5;            % 1:5
        1:6};           % 1:5

% convention:
% i step index
% j case index
% k file index

%% Script body
%% 1rename files
%matlabpool open 6

if any(itr_stp==1)
    tic
    s2_rename
    toc
end
matlabpool close

%% 2 filter images 
% max 3 worker, 1h@home 10h@psi
if any itr_stp==2
    tic
    s2_ImageFilter
    toc
end
%% 3 add images
% can use mor than 4 worker, maybe 6. 10min@home@4worker

tic
s2_ImageAdder
toc
matlabpool close

%% 4 correct for DC, OP & reference
tic 
s2_corr
toc


%% 5 make sinograms
s2_sinograms
% negative values occur in sinograms 2049 and 2050


%% 6 centering

%% 



; %1

; %2

; %3

%5

%6

%7

%8

%9

%10
