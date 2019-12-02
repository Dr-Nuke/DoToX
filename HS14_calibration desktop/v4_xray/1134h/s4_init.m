
clc             % clear commad line output
clearvars       % clear variables
clear global    % clear globals
clear all
close all       % close figures
format compact  % compact command line output



% add path so i can use the same functions
% win
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')
scriptpath='C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3/';
dir_orig ='C:/data/Tomo_HS14/02_rawdata/';
dir_work ='C:/data/Tomo_HS14/processed/';

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