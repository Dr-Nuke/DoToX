% this is the 3rd version of total,
% for the slice wise processing and intermediate saving to hard disk

clc             % clear command line output
clearvars       % clear variables
clear global    % clear globals
format compact  % compact command line output

%get name of thisscript
[~,name,~]=fileparts(mfilename('fullpath')); %get m file name

% add path so i can use the same functions
scriptpath=pwd;
% win
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')
dir_orig ='C:\data\Tomo_HS14\02_rawdata\';
dir_work ='C:\data\Tomo_HS14\processed\';

% redhat
% addpath('/home/cbolesch/HS14_calibration desktop')
% addpath('/home/cbolesch/HS14_calibration desktop/v2_1')
% addpath('/home/cbolesch/HS14_calibration desktop/v3')
% dir_orig ='/media/data/cbolesch/Tomo_HS14/02_rawdata/';
% dir_work ='/media/data/cbolesch/Tomo_HS14/processed/';

% pre- settings
t_tot=tic;
tic;

% time message string
t_string='%s scripttime: %.0f seconds,; total time %.0f minutes or %.0f hours\n';
%maximum pixel numbers (size of raw files)
imrange=[1,1200;  % x min and max number of pixel of raw images
         1,2450]; % y min and max number of pixel of raw images

impix={imrange(1,1):imrange(1,2);  %list of pxel
        imrange(2,1):imrange(2,2)};

% number of images per projection
num_proj=4;    
    
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
nangles=376; %number of angles
angles=(linspace(0,360,nangles));   % list of angles
angles=angles(1:end-1);             % remove 360�=0�

%LFT finding parameters
cam_res=22.1;           % pixel/mm, camera characteristic
r_rod=5.14;             % physical radius of the heating rod




% pixel window for dose correction calculation. needs to be outside of any
% structures on all the images

xdo_1 = 300; xdo_2=1050;
xdose=[1:xdo_1-imrange(1,1),xdo_2-imrange(1,1):imrange(1,2)-imrange(1,1)]; 
ydose=1:length(impix{2});

dir_raw={'empty/','D2O/','chcl3/','CD/'};            % original file paths
nam_raw={{'empty_'},{'D2O_','D2O_r_'},{'CHCl3_'}};       % raw file names
nam_xls={'empty.xls','D2O.xls','CHCL3.xls','DC_fake.xls'};        % xls file names
pfx_cas={'emp_','d2o_','cl3_','dc__'};      % work file name prefixes
dir_cas={'emp/','d2o/','cl3/','dc_/'};        % sub cas folders
dir_stp={'1raw/'};%,'2flt/','3add/','4cor/','5sin/'};  % step folders
nam_stp={'1raw','2blc','3flt','4add','5log','6tilt','7rec','8lft'}; 
log_nam={'log_emp','log_d2o','log_cl3','log_dc_'}; %
nam_block='block_';
nam_DC='DC';
nam_log='logfile';
nam_log_ext=strcat(nam_log,'.mat');
nam_log_ful=strcat(dir_work,nam_log_ext);


nfiles={{1:406}...          % 1:406 %how does this work?
        {1:938,939:1619}...  % {1:938,1:681}  mind special 
        {1:1621}...         % 1:1621
        {1:5}...            % 1:5
        {1:6}};           % 1:6
nfilesmax=[406,1619,1621,5,6];

ind_jj=[1,2,1,1,1]; % needed for d2o special treatment

tiltline=[600, 2000]; % lines for tilt correction

% block_size = 400; %specifies the data block size (memory limitation)
% used for low memory machines

% index convention
% i step
% j case, jj for d20 subcase
% k file number 
% itr_cas_emp=[1]; %max is 3

% step selection
% which of the steps should be performed? full scrip is [1,2,3,4,5]
itr_stp=[1,2,3,4,5,6,7,8];    % 1-image sorting, 2- block building
                        % 3-filtering
itr_stp=[];

% case selection
itr_cas=[1,2,3,4]; % 1-empty, 2-d2o, 3-cl3, 4-dc
%itr_cas=[]

% filter selection
itr_flt=[1,2,3,4]; % 1- spot 2-DC 3-Dose 4-Ring

%itr_flt=[];
% read from / write to diskfile or just calculate them:
ind_write=1; %
ind_read=1;


%% Script body

%% 1 rename and sort files
s3_rename2

%% 2 arange into blocks for paralle computing

s3_BlockForming

%% 3 apply filters
s3_FiltersDC
s3_Filters

%% 4 add images
s3_Adder

%% 5 log 
s3_divide

%% 6 tilt correction
s3_tilt  

%% 7 centereing
s3_centering

%1

%2

%3

%5

%6

%7

%8

%9

%10

t2=toc;
fprintf(t_string,name,t2-t_tot,t2/60,t2/3600)





















%