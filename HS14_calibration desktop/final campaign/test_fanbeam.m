% first load some data
 load('T')
 %%
 d=f.BoLoad(T.fnames.add{6},T);
 dat=d(173:491,:,:);
 %%
 sinogram=-log(squeeze(dat(:,640,:))); % the middle plane
sinogram=f.fraccircshift(sinogram,-(T.Cen.fitshift(2,640)-0.5))';

close all
figure(1);clf
imshow(sinogram,[])
tA.volpix=0.127; % width of one voxel in the reconstruction area
tA.detpix=0.127; % detector pixel size
tA.xPix=319; % number of x-pixels
tA.ang=linspace(0,2*pi,1357+1);
tA.ang(end)=[];
tA.orig=950/tA.detpix; %source - axis distance
tA.det=50/tA.detpix; % detector - axis distance




% create geometries and projector
vol_geom = astra_create_vol_geom(tA.xPix,tA.xPix);

% proj_geom = astra_create_proj_geom('fanflat', tA.detpix, tA.xPix, tA.ang,...
%     tA.orig, tA.det);

proj_geom = astra_create_proj_geom('fanflat', 1, tA.xPix, tA.ang,...
    tA.orig, tA.det);

%proj_id = astra_create_projector('linear', proj_geom, vol_geom);

% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);


% reconstruct
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
cfg = astra_struct('FBP_CUDA');
%cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
fbp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fbp_id);
V = astra_mex_data2d('get', recon_id);

figure(2);clf
imshow(V, []);

% garbage disposal
astra_mex_data2d('delete', sinogram_id, recon_id);
%astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', fbp_id);