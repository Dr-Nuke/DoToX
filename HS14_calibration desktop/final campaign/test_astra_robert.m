
close all

clc
load sino1.mat
sino1 = imrotate(sino1,90); % rotate 90 degrees to switch projection/pixel axes
sino1 = im2double(sino1); % make data type double
sg = sino1; % rename to sg to shorten it
%sg=nsc_sinogram(:,:,vind)
% shift rotation so that it's kind of square w.r.t. the image outline
nshift = 128; pt1 = sg(:,nshift+1:end); pt2 = sg(:,1:nshift); sg = [pt1 pt2];


ntrim = 30; sg = sg(ntrim+1:end-ntrim,:); % trim each end by ntrim pixels
[nPixels,nAngles] = size(sg); % get number of pixels and angles
stepSize = 360/nAngles; % get step size assuming equal from 0 to 360
angleListDeg = [0:stepSize:360-stepSize];
I1 = iradon(sg, angleListDeg, 'linear','Shepp-Logan',230); % I don't understand how the scaling works, the pixels don't seem to be 0.127 mm (?)
figure;
imshow(I1,[])
title('MATLAB - PARALLEL FBP')

% 0.127


I2 = ifanbeam(sg,1000/0.127,...  % sinogram, source to object center
    'FanCoverage','Cycle',... %
    'FanRotationIncrement',stepSize,...
    'FanSensorGeometry','Line',...
    'FanSensorSpacing',1,...
    'Filter','Shepp-Logan',...
    'OutputSize',230,...
    'Interpolation','linear');
figure;
imshow(I2,[])
title('MATLAB - FAN-BEAM FBP')

%%

% need to define "mySinogram" first
% add mex and tools path
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')
angleListRad = angleListDeg*pi/180; % list of projection angles
pix_per_mm = 1/0.127; % size in reconstructed image?
detectorWidth = 0.127 ; % detector pitch
reconImageSize = 230; % reconstructed image size in pixels
sourceOrigin = 950; % source to object center distance
originDetector = 50; % object center to middle detector distance
% make volume geometry defining object space
halfSize = reconImageSize/2 * 0.127;
vol_geom = astra_create_vol_geom(reconImageSize, reconImageSize, -halfSize, halfSize, -halfSize, halfSize);
% make projection geometry defining source-detector trajectories
proj_geom = astra_create_proj_geom('fanflat',...
    detectorWidth, nPixels, angleListRad, sourceOrigin, originDetector);
% projector
proj_id = astra_create_projector('strip_fanflat',proj_geom,vol_geom);
% Create the sinogram data object from sinogram data matrix
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sg');
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
I3 = astra_mex_data2d('get', rec_id); % get result
%rec = imrotate(rec,180,'crop');
figure;
imshow(I3,[])
title('ASTRA TOOLBOX - FBP')
colorbar


%%


% Create the sinogram data object from sinogram data matrix
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sg');
% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);
% Set up the parameters for a reconstruction algorithm using the GPU
%cfg = astra_struct('SIRT');
cfg = astra_struct('SIRT');
cfg.ProjectorId = proj_id;
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% perform reconstruction
nIterations = 300;
astra_mex_algorithm('iterate', alg_id, nIterations); % run algorithm
I4 = astra_mex_data2d('get', rec_id); % get result
figure;
imshow(I4,[])
title(['ASTRA TOOLBIX, SIRT (ALGEBRAIC), n=' num2str(nIterations)])
%x = mean(rec(:));
%rec(1,1)=x; rec(1,end)=x; rec(end,1)=x; rec(end,end)=x;
%v=5; % trim value
%rec(1:v,1:v)=0;
%rec(end-v+1:end,1:v)=0;
%rec(1:v,end-v+1:end)=0;
%rec(end-v+1:end,end-v+1:end)=0;
%imshow(rec, []);
%rec = astra_mex_data2d('get', rec_id); % get result
%rec = imrotate(rec,180,'crop');
    
% nIterations = 30; % iterations per step
% counter=0;
% for i=[1 2 4 8 16 32 64 128 256]
%     counter=counter+1;
%     astra_mex_algorithm('iterate', alg_id, nIterations); % run algorithm
%     rec = astra_mex_data2d('get', rec_id); % get result
%     x = mean(rec(:));
%     rec(1,1)=x; rec(1,end)=x; rec(end,1)=x; rec(end,end)=x;
%     subplot(3,3,counter) % go to subplot
%     v=5; % trim value
%     rec(1:v,1:v)=0;
%     rec(end-v+1:end,1:v)=0;
%     rec(1:v,end-v+1:end)=0;
%     rec(end-v+1:end,end-v+1:end)=0;
%     imshow(rec, []);
% end

%%

% trimming high/low vals and normlizing just for qualitative checks

I1N = RC_adjustScale(I1,1);
I2N = RC_adjustScale(I2,1);
I3N = RC_adjustScale(I3,1);
I4N = RC_adjustScale(I4,1);
I1N = RC_normalizeMatrix(I1N);
I2N = RC_normalizeMatrix(I2N);
I3N = RC_normalizeMatrix(I3N);
I4N = RC_normalizeMatrix(I4N);

figure
[y1,x1] = hist(I1N(:),100);
[y2,x2] = hist(I2N(:),100);
[y3,x3] = hist(I3N(:),100);
[y4,x4] = hist(I4N(:),100);
plot(x1,y1);
hold on
plot(x2,y2);
plot(x3,y3);
plot(x4,y4);
legend('1','2','3','4')

%%

% use normalize stuff just for visualization
figure

subplot(2,2,1)
imshow(I1N,[0 1])
title('MATLAB, PARALLEL FBP')

subplot(2,2,2)
imshow(I2N,[0 1])
title('MATLAB, FAN-BEAM FBP')

subplot(2,2,3)
imshow(I3N,[0 1])
title('ASTRA TOOLBOX, FAN-BEAM FBP')

subplot(2,2,4)
imshow(I4N,[0 1])
title(['ASTRA TOOLBIX, SIRT (ALGEBRAIC), n=' num2str(nIterations)])

%%

figure
T1 = I1N(1:100,1:100);
T2 = I2N(1:100,1:100);
T3 = I3N(1:100,1:100);
T4 = I4N(1:100,1:100);

subplot(2,2,1)
imshow(T1,[0 1])
title('MATLAB, PARALLEL FBP')

subplot(2,2,2)
imshow(T2,[0 1])
title('MATLAB, FAN-BEAM FBP')

subplot(2,2,3)
imshow(T3,[0 1])
title('ASTRA TOOLBOX, FAN-BEAM FBP')

subplot(2,2,4)
imshow(T4,[0 1])
title(['ASTRA TOOLBIX, SIRT (ALGEBRAIC), n=' num2str(nIterations)])

%%

figure
T1B = T1(51:100,1:50);
T2B = T2(51:100,1:50);
T3B = T3(51:100,1:50);
T4B = T4(51:100,1:50);

subplot(2,2,1)
imshow(T1B,[0 1])
title('MATLAB, PARALLEL FBP')

subplot(2,2,2)
imshow(T2B,[0 1])
title('MATLAB, FAN-BEAM FBP')

subplot(2,2,3)
imshow(T3B,[0 1])
title('ASTRA TOOLBOX, FAN-BEAM FBP')

subplot(2,2,4)
imshow(T4B,[0 1])
title(['ASTRA TOOLBIX, SIRT (ALGEBRAIC), n=' num2str(nIterations)])

