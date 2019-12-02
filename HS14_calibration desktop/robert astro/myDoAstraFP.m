

I = imrotate(I,180);
volume_id = astra_mex_data2d('create', '-vol', vol_geom, I);
% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP');
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id);
sinogram = sinogram';

