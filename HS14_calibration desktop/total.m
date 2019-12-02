% This script is the total cordination of the reconstruction from 
% the measurements of HS14, which featued a static setup without
% heating and flow.

%% todo

% done: put fitsread & gammafilter into function
% 2) consider zeros in late cl images , starting at cl3_1104 (eventually
% only above the plug)
% 3) cl3 OB is screwed, atm the late cl3 pictures are brighter in the right
% half than left half
% done zero and negative correction
% 5) there is a definite border between cl3 images and degraded cl3 images
% 6) gefilterte Bilder speichern
% done DC bzw scatter correction
% 8) bessere zero-behandlung als "if <=0 then put to 0.001" in f_DC
% done: line filter für cl3 projektionen ab 258
% done: mask generator in s_recon is buggy (test mit n-rec =9)
% 11 custom median filter: die eckpixel sind noch schwarz
% 12 add maxlag to xcv in s_rec
% 13 the path endpoints: the center point can move closer to the wall

%% Fragen:
% 1) centering: wenn d2o um 100 pixel verschoben wird und cle um 101, was
% dann?
% 2) wie normieren, mit den 1.0xx - farbwerten?


%%
clc
clear all
format compact


% remarks:
% one image corresponds to 20 seconds of exposure

%% pre- settings

% directories
tic
readpath='C:\data\tomo_HS14\02_rawdata\';
writepath='C:\data\tomo_HS14\processed\';
scriptpath='V:\1phd\codes\reconstruction\HS14_calibration\';

%maximum pixel numbers (size of raw files)
xmin=91; ymin=1;%consider cropping 
xmax=1200; ymax=2450; 

ind_debug=1; % set to 0 to save memory, set to 1 to have more workspace data active


yminn=[500:100:2100,2201];% in steps of 100


% threshold for median filter [if (raw-filter) > thresh_medi*raw => correct]
thresh_medi=0.02;
;
% pixel resolution
res=22.1; %pixel/mm, found in resolution.m
%%
for iii=12:12%length(yminn)-1
    ymin=yminn(iii);ymax=yminn(iii+1)-1;
    fprintf('\n running pixel rows %d to %d \n', ymin, ymax);

% image cropping for fitsread. to save resources while developping
%xpix=[xmin:xmax]; ypix=[ymin,335:355,1200:1210,2115:2140,ymax]; % 90: bad region of DC images
xpix=[xmin:xmax]; ypix=[ymin:ymax];
% angle list
angles=(linspace(0,360,376));angles=angles(1:end-1);

% pixel window for dose correction calculation. needs to be outside of any
% structures on all the images
xdo_1 = 300; xdo_2=1050;
xdose=[1:xdo_1-xmin,xdo_2-xmin:xmax-xmin]; 
ydose=[1:length(ypix)];

%% 1) empty channel
run(strcat(scriptpath,'s_empty'))


%% 2) D2O
run(strcat(scriptpath,'s_d2o'))


%% 3) CHCl3
run(strcat(scriptpath,'s_cl3'))


%% 5) DC
run(strcat(scriptpath,'s_DC'))


%%
% save('C:\data\tomo_HS14\processed\workspace.mat',...
%     'empty','empty_180','empty_OB','d2o','d2o_180','d2o_OB',...
%     'cl3','cl3_180','cl3_OB')

%% 4) corrections 
run(strcat(scriptpath,'s_correct'))


%% 5) calculations
run(strcat(scriptpath,'s_calc'))


%% 6) Reconstruction
run(strcat(scriptpath,'s_recon'))

%% 7) save files
%run(strcat(scriptpath,'s_save'))


end

