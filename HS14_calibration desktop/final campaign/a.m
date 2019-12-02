classdef a % Astra reconstruction
methods(Static)
    
function A=CreateAstraDataStruct()
    % creates the properties for the astra struct
    %create angles
    ang=linspace(0,2*pi,1357+1);
    ang(end)=[];
    A.ang=ang;  % angle list
    A.Dpix=0.127; % recon size
    A.Vpix=0.127; % recon size
    A.xPix=319; % number of x-pixels
    
end
    
function rec=FBP(sino,A)
    % FBP parallel strip
    
    vol_geom = astra_create_vol_geom(A.xPix, A.xPix);
    proj_geom = astra_create_proj_geom('parallel', 1,...
        A.xPix, A.ang);
    proj_id = astra_create_projector('strip', proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
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
    
end

function rec=FBPexpl(sino,xPix,ang)
    % FBP parallel strip
    
    vol_geom = astra_create_vol_geom(xPix, xPix);
    proj_geom = astra_create_proj_geom('parallel', 1,...
        xPix, ang);
    proj_id = astra_create_projector('strip', proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
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
    
end

function rec=FBPexplFan(sino,recImSize,ang,detPitch,src,det)
    % FBP fanbeam strip
    % sino Sinogram
    % recImSize % reconstructed image size in pixels
    % ang the angles array in radian
    % src source-axis distance
    % detPtch  detector pixel spacing
    % det stector-axis distance
    
    halfSize = recImSize/2*detPitch;
    xPix=size(sino,2);
    
    vol_geom = astra_create_vol_geom(recImSize, recImSize, ...
       -halfSize, halfSize, -halfSize, halfSize);
    
    proj_geom = astra_create_proj_geom('fanflat',...
        detPitch, xPix, ang, src, det);

%     proj_id = astra_create_projector('strip_fanflat',...
%         proj_geom, vol_geom);
    proj_id = astra_create_projector('line_fanflat',...
        proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
    rec_id = astra_mex_data2d('create', '-vol', vol_geom);
    
    % create configuration
    cfg = astra_struct('FBP_CUDA');
    %cfg.ProjectorId = proj_id;
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;
    % cfg.FilterType = 'Ram-Lak';
    
    alg_id = astra_mex_algorithm('create', cfg);
    
    astra_mex_algorithm('run', alg_id);
   
    rec = astra_mex_data2d('get', rec_id);
    astra_mex_algorithm('delete', alg_id);
    astra_mex_data2d('delete', rec_id);
    astra_mex_data2d('delete', sinogram_id);
    
    
end



function rec= SIRTexpl(sino,xPix,ang,iter)


    
    % create geometries
    vol_geom = astra_create_vol_geom(xPix, xPix);
    proj_geom = astra_create_proj_geom('parallel', 1,...
        xPix, ang);
    
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
    % reconstruct
    recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    cfg = astra_struct('SIRT_CUDA');
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;
%     cfg.option.MinConstraint = 0;
%     cfg.option.MaxConstraint = 255;
    sirt_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('iterate', sirt_id,iter);
    rec = astra_mex_data2d('get', recon_id);

    
    % garbage disposal
    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_algorithm('delete', sirt_id);
end

function rec=SARTexpl(sino,xPix,ang,iter)
    
    % create geometries
    vol_geom = astra_create_vol_geom(xPix, xPix);
    proj_geom = astra_create_proj_geom('parallel', 1,...
        xPix, ang);
    
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    % reconstruct
    recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    cfg = astra_struct('SART_CUDA');
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;
%     cfg.option.ProjectionOrder = 'custom';
%     cfg.option.ProjectionOrderList = [0:5:175 1:5:176 2:5:177 3:5:178 4:5:179];
    sart_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('iterate', sart_id, iter);
    rec = astra_mex_data2d('get', recon_id);
     
    % garbage disposal
    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_algorithm('delete', sart_id);
end

function rec=FBPexplCone(sino,volxpix,volzpix,ang,detPitch,src,det)
    % FBP conebeam strip
    % sino Sinogram
    % recImSize % reconstructed image size in pixels
    % ang the angles array in radian
    % src source-axis distance
    % detPtch  detector pixel spacing
    % det stector-axis distance
    
    halfSize = recImSize/2*detPitch;
    Volx=volxpix*detPitch;
    Voly=volxpix*detPitch;
    Volz=volzpix*detPitch;
    xPix=size(sino,2);

    vol_geom = astra_create_vol_geom(Volx,Voly,Volz);
    
    proj_geom = astra_create_proj_geom('cone',...
        detPitch,detPitch, volzpix,...
        volxpix, ang, src, det);
    
    rec_id = astra_mex_data2d('create', '-vol', vol_geom);
    
    proj_id=astra.data3d.create('-proj3d',proj_geom,P)
    rec_id=astra.data3d.create('-vol',vol_geom)
    
    
    proj_id = astra_create_projector('strip_fanflat',...
        proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
    
    
    % create configuration
    cfg = astra_struct('FBP_CUDA');
    %cfg.ProjectorId = proj_id;
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;
    %cfg.FilterType = 'Ram-Lak';
    
    alg_id = astra_mex_algorithm('create', cfg);
    
    astra_mex_algorithm('run', alg_id);
   
    rec = astra_mex_data2d('get', rec_id);
    astra_mex_algorithm('delete', alg_id);
    astra_mex_data2d('delete', rec_id);
    astra_mex_data2d('delete', sinogram_id);
    
end




end
end

