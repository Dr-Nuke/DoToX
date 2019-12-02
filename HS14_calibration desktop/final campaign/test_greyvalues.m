
plane=600;
mysino=squeeze(d(:,plane,:));

% rec=a.FBPexplFan(...
%     mysino',...
%     recsize,ang,detPitch,src,det);
%%
halfSize = recsize/2*detPitch;
xPix=size(mysino,1);

vol_geom = astra_create_vol_geom(recsize, recsize, ...
    -halfSize, halfSize, -halfSize, halfSize);

proj_geom = astra_create_proj_geom('fanflat',...
    detPitch, xPix, ang, src, det);


proj_id = astra_create_projector('line_fanflat',...
    proj_geom, vol_geom);
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, mysino');

rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% create configuration
cfg = astra_struct('FBP_CUDA');
%cfg.ProjectorId = proj_id;
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
% cfg.FilterType = 'Ram-Lak';

alg_id = astra_mex_algorithm('create', cfg);

astra_mex_algorithm('run', alg_id);

recfan = astra_mex_data2d('get', rec_id);
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
%%
    vol_geom = astra_create_vol_geom(xPix, xPix);
    proj_geom = astra_create_proj_geom('parallel', 1,...
        xPix, ang);
    proj_id = astra_create_projector('line', proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, mysino');
    
    rec_id = astra_mex_data2d('create', '-vol', vol_geom,0);
    
    % create configuration
    cfg = astra_struct('FBP_CUDA');
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;
    cfg.FilterType = 'Ram-Lak';
    
    alg_id = astra_mex_algorithm('create', cfg);
    
    astra_mex_algorithm('run', alg_id);
   
    rec = astra_mex_data2d('get', rec_id);
    astra_mex_algorithm('delete', alg_id);
    astra_mex_data2d('delete', rec_id);
    astra_mex_data2d('delete', sinogram_id);
    
    %%

figure(1);clf;
imshow(mysino,[])

figure(2);clf;
imshow(rec,[])
title('parallel')

figure(3);clf;
histogram(rec,100)
title('parallel')


figure(4);clf;
imshow(recfan,[])
title('fan')
figure(5);clf;
histogram(recfan,100)
title('fan')


alfan=8.3e-6;
alpar=0.00734;
