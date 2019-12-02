
%%

% minimal working example for phantom-forwardprojection-reconstruction
% via astra fanbeam

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


%% parameters

sial=0.9; % aluminum attenuation coefficient 1/cm
d=10;% phantom area width in mm
n=100; % Phantom area width in pixel [for 0.1 mm pixels]
pitch_f=d/n; % phantom pixel pitch in mm
ndet=200;       % number of detector pixels
detpitch=0.1;   % detector pixel spacing in mm
recpitch=0.1;   % reconstruction pixel size in mm
drec=10;        % reconstruction area width in mm
nrec=round(drec/recpitch); % reconstruction area width in pixel
nang= 1000; % number of projection angles
src=950;   % source-origin-distance in mm
det=50;    % detector-origin-distance in mm
% phantom
V = zeros(n);
gv=sial/10; % alu attenuation in 1/mm
gv=1;       % no^^ we stick to minimal working example values
V(25:75,25:75)=gv; % set some phantom part to unity
angtest=linspace(0,2*pi,nang+1); % make the angle list in 2 steps
angtest(end)=[];

%% create geometries

proj_geom = astra_create_proj_geom('fanflat', detpitch ,ndet, angtest,...
    src,det);
vol_geom = astra_create_vol_geom(n,n,-d/2,d/2,-d/2,d/2);

%% store volume

volume_id = astra_mex_data2d('create', '-vol', vol_geom, V);

%% create forward projection


proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
cfg.ProjectorId = proj_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id); % this appears to have units of 
astra_mex_data2d('delete', sinogram_id, volume_id);
astra_mex_algorithm('delete', fp_id);

%% plots

fi=1;
figure(fi);clf;fi=fi+1; clf; subplot(1,2,1);
imshow(V,[]);set(gca,'YDir','normal')
title('phantom')
colorbar

fi1=fi;
figure(fi);clf;fi=fi+1;
subplot(1,2,1)
imshow(sinogram, []);set(gca,'YDir','normal')
colorbar
title('$\int \Sigma ds$ (sinogram)','Interpreter','latex')
xlabel(sprintf('%f',max(sinogram(1,:))))
sino2=exp(-sinogram);
subplot(1,2,2)
imshow(sino2, []);set(gca,'YDir','normal')
colorbar
title('$e^{-\int \Sigma ds}$ (measurement)','Interpreter','latex')


%% recon


%sinogram = sinogram/pitch_f; % RA: this line converts sinogram values from 1/mm to 1/pixel  ????????????????????????????????

xPix=size(sinogram,2); % sinogram width; here equal to the detector width
vol_geom2 = astra_create_vol_geom(nrec, nrec, ...
    -drec/2, drec/2, -drec/2, drec/2);
%proj_geom  = astra_create_proj_geom('fanflat', detpitch ,ndet, angtest, src, det); % RA: copied here in comment for comparison
 proj_geom2 = astra_create_proj_geom('fanflat', detpitch, xPix, angtest, src, det);
%sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);  % RA: copied here in comment for comparison
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom2, sinogram);
rec_id2 = astra_mex_data2d('create', '-vol', vol_geom2);

proj_id2 = astra_create_projector('line_fanflat', proj_geom2, vol_geom2);

% create configuration
cfg2 = astra_struct('FBP_CUDA');
cfg2.ProjectorId = proj_id;
cfg2.ReconstructionDataId = rec_id2;
cfg2.ProjectionDataId = sinogram_id2;
%cfg.FilterType = 'Ram-Lak';
alg_id2 = astra_mex_algorithm('create', cfg2);
astra_mex_algorithm('run', alg_id2);
% get the reconstruction
rec = astra_mex_data2d('get', rec_id2);%/pitch_f % RA: converted from 1/pixel to 1/mm ?????????????????????????????????????
astra_mex_algorithm('delete', alg_id2);
astra_mex_data2d('delete', rec_id2);
astra_mex_data2d('delete', sinogram_id2);

figure(1); subplot(1,2,2);
imshow(rec,[]);set(gca,'YDir','normal')
colorbar
title('recon')