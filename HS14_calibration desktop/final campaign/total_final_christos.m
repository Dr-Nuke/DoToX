%%

% this is the script that evaluates the data from 
% the final measurement campaign in Descember 2017

% data is being held on disk, parameters are saved in struct T (Total)
% T total struct, frequently used parameters go here
% T.d data (raw) related properties
% T.dose contains dose corection related stuff

%% initialization
clc             % clear commad line output
clear global    % clear globals
clear           % clear workspace
close all       % close figures
format compact  % compact command line output
I1=f_C; %I1 for Instance 1
I1.T.sys.t_tot=tic;
tic;
fprintf('running the tomo processing script...\n')

% .d stands for data, raw data, and related properties

I1.T.d.DataPath='E:\20171213 final campaign\'; % dir of raw data

%T.d.DataPath='C:\data\final campaign\';
I1.T.proofs=1; % flag to activate (or deactivate) prooving routines & user interactivities
I1.T.d.header=2048; % the file header to be skipped at read in
I1.T.d.imsize=[640 1024]; % the frame size in pixels
I1.T.d.GreyBytes=2; % the color depth in bytes, usually 2 for VIVA .seq videos
I1.T.d.ncas=9; %number of cases/ Betriebspunkte
I1.T.d.nrep=10; % number of repetitions for each case


% generate list of raw data files
I1.T.d.fileflag=('*tom_*.seq'); % files with this flag will enter the file list
I1.GenerateFileList();    % finally read in the raw files list

