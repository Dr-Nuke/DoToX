% test centering by upsampling

% prepare data
close all
clc
clear

load('T.mat')
%%
cas=6;
rep=1;
T.cas=6;
T.rep=1;

sinsh=T.Rec.SinoRotShift;
ang=deg2rad(T.Rec.angles);
recsize=T.Rec.recsize;
detPitch=T.Rec.detPitch;
src=T.Rec.src;
det=T.Rec.det;
niter=T.Cen.niter;
BS=T.Raw.BS;
xmin=T.Cen.xmin;
xmax=T.Cen.xmax;
MaskThresh=T.Cen.MaskThresh;
recImSize=recsize;

fname=sprintf('%s%02d_%02d.mat','1_corr_',T.cas,T.rep);
fpath=sprintf('%s%s',T.d.DataPath,fname);
d=f.BoLoad(fpath,T);

% crop z
d=d(:,:,T.q360.startframe:(...
    T.q360.startframe+T.q360.nFrames(1,1)-1));

% image correction: pixels & filter
corr=f.ImCorrFilt(mean(d,3),T);
% apply image correction
d=f.ApplyImCorr(d,corr,T);

% crop x & y
d=d(T.Cen.range,:,:);

% quadrant-rotate, fliplr, log
d=flip(-log(circshift(d,T.Cen.quadrot,3)),3);

% find the rough centering for each frame ~+-2pix
[T.Cen.fit{cas,rep},T.Cen.fitshift(cas,rep,:),T.Cen.centshift(cas,rep,:)]=...
    f.FindCentering(T,d);

baseshift=squeeze(T.Cen.fitshift(cas,rep,:));
% if all([cas,rep]==[1,1])
T.Cen.MMask(:,:,cas,rep)=f.MakeTomoMask('par',f.fraccircshift(...
    squeeze(d(:,T.Cen.MaskPlane,:)),...
    -baseshift(T.Cen.MaskPlane)),recsize,ang,detPitch,src,det,...
    MaskThresh,0.000005,[150 150]);
f.CheckMMask(T.Cen.MMask(:,:,cas,rep),cas,rep)
% end
MMask=T.Cen.MMask(:,:,cas,rep);
%%
sf=10;
niter=(sf)*4+1;
smin=-2*sf;
smax=2*sf;
shift=linspace(smin,smax,niter);
planes=[100,350,512,1000];
planes=10:1000;

q=nan(niter,length(planes));
xlist=nan(niter,length(planes));
eshift=nan(1,length(planes));


xPix=recsize;
%rec=zeros(recsize,recsize(niter,length(planes)))
for p=1:length(planes);
    plane=planes(p);
    disp(plane);
    sino=squeeze(d(:,plane,:));
    % upsample
    us=size(sino).*[sf,1];
    %sino2=imresize(sino,us,'bilinear');
    sino2=imresize(sino,us,'bilinear');
    if plane==10; % fist case
        for s=1:niter
            
            sino3=circshift(sino2,shift(s),1);
            vol_geom = astra_create_vol_geom(xPix, xPix);
            proj_geom = astra_create_proj_geom('parallel', 1/sf,...
                us(1), ang);
            proj_id = astra_create_projector('strip', proj_geom, vol_geom);
            sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino3');
            
            rec_id = astra_mex_data2d('create', '-vol', vol_geom,0);
            
            % create configuration
            cfg = astra_struct('FBP_CUDA');
            cfg.ReconstructionDataId = rec_id;
            cfg.ProjectionDataId = sinogram_id;
            cfg.FilterType = 'Ram-Lak';
            
            alg_id = astra_mex_algorithm('create', cfg);
            
            astra_mex_algorithm('run', alg_id);
            
            rec2 = astra_mex_data2d('get', rec_id);
            astra_mex_algorithm('delete', alg_id);
            astra_mex_data2d('delete', rec_id);
            astra_mex_data2d('delete', sinogram_id);
            [q(s,p),~]=qc.VarianceQuality(rec2,MMask);
            %         figure();
            %         imshow(rec2,[])
            %         title(num2str(shift(s)))
        end
    else
        for s=(vind(p-1)-2):(vind(p-1)+2)
            
            sino3=circshift(sino2,shift(s),1);
            vol_geom = astra_create_vol_geom(xPix, xPix);
            proj_geom = astra_create_proj_geom('parallel', 1/sf,...
                us(1), ang);
            proj_id = astra_create_projector('strip', proj_geom, vol_geom);
            sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino3');
            
            rec_id = astra_mex_data2d('create', '-vol', vol_geom,0);
            
            % create configuration
            cfg = astra_struct('FBP_CUDA');
            cfg.ReconstructionDataId = rec_id;
            cfg.ProjectionDataId = sinogram_id;
            cfg.FilterType = 'Ram-Lak';
            
            alg_id = astra_mex_algorithm('create', cfg);
            
            astra_mex_algorithm('run', alg_id);
            
            rec2 = astra_mex_data2d('get', rec_id);
            astra_mex_algorithm('delete', alg_id);
            astra_mex_data2d('delete', rec_id);
            astra_mex_data2d('delete', sinogram_id);
            [q(s,p),~]=qc.VarianceQuality(rec2,MMask);
            %         figure();
            %         imshow(rec2,[])
            %         title(num2str(shift(s)))
        end
    end
    
    [vmax(p),vind(p)]=max(q(:,p));
    sino3=circshift(sino2,shift(vind(p)),1);
    halfSize = recImSize/2*detPitch;
    
    vol_geom = astra_create_vol_geom(recImSize, recImSize, ...
        -halfSize, halfSize, -halfSize, halfSize);
    
    proj_geom = astra_create_proj_geom('fanflat',...
        detPitch/sf, us(1), ang, src, det);
    
    proj_id = astra_create_projector('strip_fanflat',...
        proj_geom, vol_geom);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino3');
    
    rec_id = astra_mex_data2d('create', '-vol', vol_geom);
    
    % create configuration
    cfg = astra_struct('FBP_CUDA');
    %cfg.ProjectorId = proj_id;
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;
    % cfg.FilterType = 'Ram-Lak';
    
    alg_id = astra_mex_algorithm('create', cfg);
    
    astra_mex_algorithm('run', alg_id);
    
    rec(:,:,p) = astra_mex_data2d('get', rec_id);
    astra_mex_algorithm('delete', alg_id);
    astra_mex_data2d('delete', rec_id);
    astra_mex_data2d('delete', sinogram_id);
    
    fid=figure(4853);clf
    fid.Position=[100 162 400 800];
    subplot('Position',[0.05,0.5,0.9,0.45]);
    imshow(rec(:,:,p),[]);
    title(sprintf('c%d r%d p%03d',cas,rep,plane));
    subplot('Position',[0.1,0.1,0.8,0.3]);
    pind=~isnan(q(:,p));
    plot(shift(pind),q(pind,p),'xr','Displayname','upsampling')
    hold on
    
    plot(shift(vind(p)),q(vind(p),p),'+b','Displayname',...
        sprintf('%.3f',q(vind(p),p)))
    title(f.FigTName('Another mask',0,cas,rep))
    title(sprintf('shift finding, %d iterations',length(pind)))
    grid on
    xlim([smin,smax]);
    legend('Location','South')
    
    
    fname=f.FigFileName3('upsampling',0,cas,rep,plane);
    f.f_BoFig3PDF(fid,'reconcheck\',fname)
    
    
    
end

































