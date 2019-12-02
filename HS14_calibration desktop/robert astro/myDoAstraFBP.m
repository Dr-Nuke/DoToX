
% Create the sinogram data object from sinogram data matrix
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, mySinogram');
% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);
% Set up the parameters for a reconstruction algorithm using the CPU
cfg = astra_struct('FBP');
cfg.ProjectorId = proj_id;
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% perform reconstruction
astra_mex_algorithm('run', alg_id); % run algorithm
rec = astra_mex_data2d('get', rec_id); % get result
rec = imrotate(rec,180,'crop');