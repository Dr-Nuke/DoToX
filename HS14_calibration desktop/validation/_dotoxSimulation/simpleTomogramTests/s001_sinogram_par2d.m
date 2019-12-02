%close all
clear
clc

% add mex and tools path for either personal laptop or ETH PC
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end


% -----------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
% 
% Copyright: 2010-2016, iMinds-Vision Lab, University of Antwerp
%            2014-2016, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://www.astra-toolbox.com/
% -----------------------------------------------------------------------

phantomSizePix = 100; % phantom size N, for NxN pixel phantom
phantomPixelSize = 0.5;

% Create a basic 100x100 square volume geometry
vol_geom = astra_create_vol_geom(phantomSizePix, phantomSizePix, ...
    -phantomSizePix/2, phantomSizePix/2, -phantomSizePix/2, phantomSizePix/2);
% Create a parallel beam geometry with 360 angles between 0 and 2*pi,
% and 100 detector pixels of width 1.
proj_geom = astra_create_proj_geom('parallel', 1.0, 100, linspace2(0,2*pi,360));
P = zeros(phantomSizePix,phantomSizePix);
P(25:75,25:75) = 1; % make middle square phantom have values of 1
% Create a sinogram using the GPU.
% Note that the first time the GPU is accessed, there may be a delay
% of up to 10 seconds for initialization.
[sinogram_id, sinogram] = astra_create_sino_gpu(P, proj_geom, vol_geom);

clf
subplot(2,2,1); imshow(P, []);
subplot(2,2,2); imshow(sinogram, []);
subplot(2,2,3); hist(sinogram(:),100); xlim([2 inf])
astra_mex_data2d('delete', sinogram_id); % Free memory
