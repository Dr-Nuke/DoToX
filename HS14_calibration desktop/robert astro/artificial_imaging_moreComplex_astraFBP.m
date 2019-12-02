
close all
clear
clc


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% NOTE!! SUB-CHANNEL VOID FRACTION NOT UPDATED TO ACCOUNT FOR COMPLEX GEOMETRY !!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% add mex and tools path for Astra toolbox and prepare data for FP
% ETH PC
addpath('C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\mex')
addpath('C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\tools')
% personal laptop
%addpath('C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\mex')
%addpath('C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\tools')

% define parameters

% for DD fast neutrons (2.8 MeV)
% Zirconium density 5.68 g/cm3, micro XS for DD = 4 barn, mass 91.2...
% gives 0.15/cm for Zr  (from data in above line)
%       0.19/cm for H2O (from thesis)
% for photons, can consider
% X-ray ~100 keV equivalent (flux from spekcalc)
% Ir-192 310 keV (141% emission)
% Ir-192 475 keV ( 52% emission)
% Ir-192 604 keV ( 18% emission)
% Cs-137 662 keV ( 80% emission)
% Co-60 1250 keV (200% emission)
E = 1.250; % Energy in MeV;
% pixel size for producing fine ideal images before compressing them
% ...must evenly divide into pitch (e.g. 0.1, 0.05, 0.025 ...)
pixelSize1 = 0.1; % [mm]
pixelSize2 = 0.25; % pixel size of compressed image used for forward projection [mm]
pixelSize3 = 0.5; %  pixel size of the reconstructed image [mm]
detectorResolution = 0.25; % detector resolution [mm]
% get attenuation for two materials at energy E [MeV]
[attnFlow,attnPin]= my_getAttenuation(E); 
attn3 = attnPin*pixelSize2/10; % attenuation of box materials is attn3
%attnPin = attnFlow; % this rewrites the pins to be water also (approx = plastic)
%attnPin = 0.15; attnFlow = 0.19; attn3 = 0.15*pixelSize2/10;  % rewrite to use D-D values
%attn3 = 0.19*pixelSize2/10;
%attnPin = 1.0; attnFlow = 0.16; attn3 = 1.0*pixelSize2/10;  % rewrite to use rough X-ray values
attnPin = 0.16; attnFlow = 0.16; attn3 = 1.0*pixelSize2/10; % rewrite to use rough X-ray values with plastic pins

% these are defined in my_makeUnitCell:
% pitch = 12.9; % mm (rounded from 12.87)
% pinOuter = 9.8; % mm
% pinInner = 8.5; % mm , from 0.65 mm cladding 
projectionStep = 1; % angle between one projection and the next
nDetectors = round(200/detectorResolution); % object diagonal a bit less than 200 mm for now
sourceObject = 1800; % mm distance source to object center
objectDetector = 200; % mm distance object center to detectors
angleSpacing = detectorResolution/(sourceObject+objectDetector)*180/pi;
angleList = [0:projectionStep:360-projectionStep]*pi/180; % list of projection angles
% this prepares the object images Iempty, Ivaried, and Iflat with pixel size pixelSize2
myMakeObjectImages2;


I=Iempty;  my_addFeatures; Iempty =I;
I=Ivaried; my_addFeatures; Ivaried=I;
I=Ifilled; my_addFeatures; Ifilled=I;

x=zeros(632,632);
x(50:532+49,50:532+49)=Ivaried;
figure;imshow(x,[]);figure;

%%

% do subchannel averaging for the ideal images
IfilledSC = imresize(Ifilled,[10 10],'box');
IemptySC = imresize(Iempty,[10 10],'box');
IvariedSC = imresize(Ivaried,[10 10],'box');
% void fraction as (full-varied)/(full-empty) calculated directly from raw image data
VF_from_original_images = (IfilledSC - IvariedSC) ./ (IfilledSC - IemptySC);

reconImageSize1 = size(Iempty,1); % take attenuation map image side length
sourceObjectPixels1 = round(sourceObject/pixelSize2); % object center to source in pixels
objectDetectorPixels1 = round(objectDetector/pixelSize2); % object center to middle detector distance
% make volume geometry defining object space based on attenuation map size
vol_geom = astra_create_vol_geom(reconImageSize1, reconImageSize1);
% make projection geometry defining source-detector trajectories
proj_geom = astra_create_proj_geom('fanflat',...
    detectorResolution/pixelSize2, nDetectors, angleList, sourceObjectPixels1, objectDetectorPixels1);
% projector
proj_id = astra_create_projector('strip_fanflat',proj_geom,vol_geom);

I = Iempty; myDoAstraFP; SG_empty = sinogram;
I = Ivaried; myDoAstraFP; SG_varied = sinogram;
I = Ifilled; myDoAstraFP; SG_filled = sinogram;

% redo some data for FBP
reconImageSize2 = round(reconImageSize1*pixelSize2/pixelSize3) ; % reconstructed image size
sourceObjectPixels2 = sourceObjectPixels1*pixelSize2/pixelSize3;
objectDetectorPixels2 = objectDetectorPixels1*pixelSize2/pixelSize3; % object center to middle detector distance
% make volume geometry defining object space
vol_geom = astra_create_vol_geom(reconImageSize2, reconImageSize2);
% make projection geometry defining source-detector trajectories
proj_geom = astra_create_proj_geom('fanflat',...
    detectorResolution/pixelSize3, nDetectors, angleList, sourceObjectPixels2, objectDetectorPixels2);
% projector
proj_id = astra_create_projector('strip_fanflat',proj_geom,vol_geom);

mySinogram = SG_empty;  myDoAstraFBP; RC_empty  = rec;
mySinogram = SG_varied; myDoAstraFBP; RC_varied = rec;
mySinogram = SG_filled; myDoAstraFBP; RC_filled = rec;

RCfilledSC = imresize(RC_filled, [10 10],'box');
RCemptySC  = imresize(RC_empty,  [10 10],'box');
RCvariedSC = imresize(RC_varied, [10 10],'box');

% void fraction as (full-varied)/(full-empty) 
% calculated directly from reconstructed image data
VF_from_reconstructed_images = (RCfilledSC - RCvariedSC) ./ (RCfilledSC - RCemptySC);
VF_error_percent = [(VF_from_original_images - VF_from_reconstructed_images)...
    ./VF_from_original_images ] *100;
VF_error_difference = [(VF_from_original_images - VF_from_reconstructed_images)] *100;
roundedErrorPercentIdeal = round(VF_error_percent*100)/100;
roundedErrorDifference = round(VF_error_difference*100)/100;
roundedErrorPercentIdeal(noPinMap==1) = 0;
roundedErrorPercentDifference(noPinMap==1) = 0;
VF_error_percent(noPinMap==1) = 0;
VF_error_difference(noPinMap==1) = 0;
VF_expected = (1-densityMap');

%%

%npsVec = [500 1000 2000 4000 8000 16000 32000 64000 128000]'
%npsVec = [100 200 400 700 1100 1600 2200 2900 3700 4700]'
%npsVec = [500 1000 2000 4000 8000]'
%npsVec = [16000 32000 64000 128000 256000]'
npsVec = [500 1000 2000 4000 8000 16000 32000 64000 128000 256000]'
%npsVec = 1000
plotCounter = 0;
ntests = 2 % number of tests to run for given nps value

for npsIndex = 1:length(npsVec);
    nps = npsVec(npsIndex)

plotCounter = plotCounter + 1;
    
resultsAll=[];

for n=1:ntests

    SG_varied_noisy = SG_varied; SG_varied_noisy(:,:) = 0;
    x = nps*exp(-SG_varied);
    SG_varied_noisy = -log( normrnd(x,sqrt(x)) / nps);
    mySinogram = SG_varied_noisy; myDoAstraFBP; RC_varied_noisy = rec;
    RC_variedSC_noisy = imresize(RC_varied_noisy,[10 10],'box');
    % void fraction as (full-varied)/(full-empty)
    % calculated directly from reconstructed image data
    VF_from_reconstructed_images_noisy = (RCfilledSC - RC_variedSC_noisy)./...
        (RCfilledSC - RCemptySC);
    VF_error_percent_noisy = [(VF_from_original_images - VF_from_reconstructed_images_noisy)...
        ./VF_from_original_images ] *100;
    VF_error_difference_noisy = [(VF_from_original_images - VF_from_reconstructed_images_noisy)] *100;
    VF_error_percent_noisy(noPinMap==1) = 0;
    VF_error_difference_noisy(noPinMap==1) = 0;
    disp(['iteration ' num2str(n) ' complete'])
    resultsAll(:,:,n) = VF_error_percent_noisy;
    resultsAll2(:,:,n) = VF_error_difference_noisy;
        
end

%%
roundedErrorPercentNoisy = mean(resultsAll,3);
roundedErrorDifferenceNoisy = mean(resultsAll2,3);
roundedErrorPercentNoisy = round(roundedErrorPercentNoisy*100)/100;
roundedErrorDifferenceNoisy = round(roundedErrorDifferenceNoisy*100)/100;

[rows,cols,layers] = size(resultsAll);
standard_deviation_of_error1 = zeros(rows,cols);
standard_deviation_of_error2 = zeros(rows,cols);
for i=1:rows
    for j=1:cols
        standard_deviation_of_error1(i,j) = std(resultsAll(i,j,:));
        standard_deviation_of_error2(i,j) = std(resultsAll2(i,j,:));
    end
end

resultsOuter(npsIndex,1) = nps;
a = standard_deviation_of_error1(:);
resultsOuter(npsIndex,2) = min(a(a>0));
resultsOuter(npsIndex,3) = mean(a(a>0));
resultsOuter(npsIndex,4) = max(a(a>0));
a = standard_deviation_of_error2(:);
resultsOuter(npsIndex,5) = min(a(a>0));
resultsOuter(npsIndex,6) = mean(a(a>0));
resultsOuter(npsIndex,7) = max(a(a>0));

subplot(2,5,plotCounter);
%I=RC_varied_noisy-RC_empty;imshow(I',[min(I(:)) max(I(:))]);
%I=RC_varied_noisy-RC_empty;imshow(I',[-1e-3 5e-3]); % gamma
%I=RC_varied_noisy-RC_empty;imshow(I',[-1e-3 12e-3]); % DD
I=RC_varied_noisy-RC_empty;imshow(I',[-0.2e-3 0.8e-3]); % xray
title(['nps = ' num2str(nps)])

end

%%

% figure;
% x = resultsOuter(:,1);
% a = resultsOuter(:,2); % min of std of error (relative)
% b = resultsOuter(:,3); % mean
% c = resultsOuter(:,4); % max
% plot(x,a,'r','LineWidth',3);
% hold on
% plot(x,b,'g','LineWidth',3);
% plot(x,c,'b','LineWidth',3);
% grid minor
% xlabel('flat field counts');
% ylabel('std dev of error (relative)');
% legend('minimum','average','maximum')
% 
% figure;
% x = resultsOuter(:,1);
% a = resultsOuter(:,5); % min of std of error (absolute)
% b = resultsOuter(:,6); % mean
% c = resultsOuter(:,7); % max
% plot(x,a,'r','LineWidth',3);
% hold on
% plot(x,b,'g','LineWidth',3);
% plot(x,c,'b','LineWidth',3);
% grid minor
% xlabel('flat field counts');
% ylabel('standard deviation of error [%]');
% legend('minimum sub-channel','average over all sub-channels','maximum sub-channel')

% figure
% subplot(1,3,1)
% scatter(VF_expected(:)*100,VF_error_percent(:),'Filled');
% xlim([80 100])
% xlabel('Void fraction [%]')
% ylabel('(Calculated-Expected)/Expected [%]')
% title('No-noise void fraction error')
% grid minor
% subplot(1,3,2)
% scatter(VF_expected(:)*100,roundedErrorPercentNoisy(:),'Filled');
% xlim([80 100])
% xlabel('Void fraction [%]')
% ylabel('(Calculated-Expected)/Expected [%]')
% title('With-noise void fraction mean error')
% grid minor
% subplot(1,3,3);
% scatter(VF_expected(:)*100,standard_deviation_of_error1(:),'Filled');
% xlim([80 100])
% xlabel('Void fraction [%]')
% ylabel('std dev of (Calculated-Expected)/Expected [%]')
% title('With-noise void fraction standard deviation of error')
% grid minor

% figure
% subplot(1,3,1)
% scatter(VF_expected(:)*100,VF_error_difference(:),'Filled');
% xlim([80 100])
% xlabel('Void fraction [%]')
% ylabel('Error [%]')
% title('No noise added')
% grid minor
% subplot(1,3,2)
% scatter(VF_expected(:)*100,roundedErrorDifferenceNoisy(:),'Filled');
% xlim([80 100])
% xlabel('Void fraction [%]')
% ylabel('Error [%]')
% title('With noise - average over all iterations')
% grid minor
% subplot(1,3,3);
% scatter(VF_expected(:)*100,standard_deviation_of_error2(:),'Filled');
% xlim([80 100])
% title('Standard deviation of noisy cases')
% xlabel('Void fraction [%]')
% ylabel('standard deviation of error [%]')
% grid minor

%%

% figure;
% subplot(2,3,1);
% I=Iempty;imshow(I,[min(I(:)),max(I(:))]);title('ideal empty')
% subplot(2,3,2);
% I=Ivaried;imshow(I',[min(I(:)),max(I(:))]);title('ideal varied')
% subplot(2,3,3);
% I=Ifilled;imshow(I,[min(I(:)),max(I(:))]);title('ideal filled')
% subplot(2,3,5);
% I=Ivaried-Iempty;imshow(I',[min(I(:)),max(I(:))]);title('ideal (varied-empty)')
% subplot(2,3,6);
% I=Ifilled-Iempty;imshow(I,[min(I(:)),max(I(:))]);title('ideal (filled-empty)')
% 
% figure;
% subplot(2,3,1);
% I=RC_varied;imshow(I,[min(I(:)) max(I(:))]);title('FBP from ideal varied')
% subplot(2,3,2);
% I=RC_varied-RC_empty;imshow(I',[min(I(:)) max(I(:))]);title('FBP from ideal (varied-empty)')
% subplot(2,3,3)
% I=RC_filled-RC_empty;imshow(I,[min(I(:)) max(I(:))]);title('FBP from ideal (filled-empty)')
% subplot(2,3,4)
% I=RC_varied_noisy;imshow(I',[min(I(:)) max(I(:))]);title('FBP with noise varied')
% subplot(2,3,5)
% I=RC_varied_noisy-RC_empty;imshow(I',[min(I(:)) max(I(:))]);title('FBP (with noise varied)-(empty)')

%figure;
%I=RC_varied_noisy-RC_empty;imshow(I',[min(I(:)) max(I(:))]);title('FBP (with noise varied)-(empty)')

