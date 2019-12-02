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

phantomSizePix = 200; % phantom matrix size N, for NxN pixel phantom
phantomObjectPix = 100; % size of square in phantom in pixels
phantomPixelSize = 1; % size of each phantom [mm]
detectorPixelSize = 1; % detector spacing [mm]
attenuation = 1; % attenuation value inserted into phantom matrix

angleList = linspace2(0,2*pi,360);  % list of projection angles
% Create a square volume geometry
C = round(phantomSizePix/2*phantomPixelSize); % half size of phantom
vol_geom = astra_create_vol_geom(phantomSizePix, phantomSizePix,-C,C,-C,C);
% Create a parallel beam geometry
proj_geom = astra_create_proj_geom('parallel', detectorPixelSize, phantomSizePix/phantomPixelSize, angleList );
P = zeros(phantomSizePix,phantomSizePix); % empty phantom matrix
% make approx. middle 10x10 pixels equal to attenuation
C1 = round(phantomSizePix/2-phantomObjectPix/2);
C2 = round(phantomSizePix/2+phantomObjectPix/2-1);
P(C1:C2,C1:C2) = attenuation; % make 5x5 pimiddle square phantom have values of 1
% Create a sinogram using the GPU.
% Note that the first time the GPU is accessed, there may be a delay
% of up to 10 seconds for initialization.
[sinogram_id, sinogram] = astra_create_sino_gpu(P, proj_geom, vol_geom);
astra_mex_data2d('delete', sinogram_id); % Free memory

% We now re-create the sinogram data object as we would do when loading
% an external sinogram
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);
% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);
% Set up the parameters for a reconstruction algorithm using the GPU
cfg = astra_struct('FBP_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
alg_id = astra_mex_algorithm('create', cfg);
% Run the algorithm
astra_mex_algorithm('run', alg_id, 150);
% Get the result
rec = astra_mex_data2d('get', rec_id);

% Clean up. (algorithm object = GPU, data objects = RAM)
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);

% show phantom, forward projected sinogram, histogram of phantom,...
% image reconstruction, histogram of reconstruction
clf
subplot(2,3,1); imshow(P, []); title('phantom');
subplot(2,3,4); hist(P(:),100); title('phantom');
subplot(2,3,2); imshow(sinogram, []); title('sinogram');
subplot(2,3,5); hist(sinogram(:),100); xlim([2 inf]); title('sinogram');
subplot(2,3,3); imshow(rec, []); title('reconstruction')
subplot(2,3,6); hist(sinogram(:),100); xlim([2 inf]); title('reconstruction');

disp(['sinogram size = ' num2str(size(sinogram,1)) ' x ' num2str(size(sinogram,2))])
