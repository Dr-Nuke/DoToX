  
function tomoNoisy = A02_doReconstruction(...
...
thisObjectImage,...
sensorResolution,...
pixelSizeCoarse,...
sensorPositions,...
projectionAngleListDeg,...
Y_nps,...
efficiency)

%% prepare stuff
% add mex and tools path for either personal laptop or ETH PC
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end
% make volume geometry defining object space
vol_geom = astra_create_vol_geom(size(thisObjectImage));
% make projection geometry defining source-detector trajectories
%proj_geom = astra_create_proj_geom('parallel',...
%    detectorWidth, nDetectors, projectionAngleListRad);
proj_geom = astra_create_proj_geom('parallel',...
    sensorResolution/pixelSizeCoarse, sensorPositions, projectionAngleListDeg*pi/180);
% projector
proj_id = astra_create_projector('line',proj_geom, vol_geom);
volume_id = astra_mex_data2d('create', '-vol', vol_geom, thisObjectImage);
% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);

%% get forward projection
astra_mex_algorithm('run', fp_id);
disp('Forward projection complete.')
sinogramClean = astra_mex_data2d('get', sinogram_id);

%% add noise to sinogram according to expected value of flat field pixels
sinogramNoisy = zeros(size(sinogramClean));
EV = [Y_nps*efficiency*exp(-sinogramClean)];

sinogramNoisy = -log(normrnd(EV,sqrt(EV))/Y_nps);

%% get noisy tomogram
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogramNoisy);
cfg = astra_struct('FBP_CUDA');% myAlg
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
fbp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fbp_id);
tomoNoisy = astra_mex_data2d('get', recon_id);
disp('Filtered back projection complete.');


%% garbage disposal
astra_mex_data2d('delete', sinogram_id, recon_id);
astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', fbp_id);

end