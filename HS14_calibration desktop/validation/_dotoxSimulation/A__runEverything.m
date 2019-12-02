

close all
clear
clc

% prepare astra toolbox stuff
% add mex and tools path for either personal laptop or ETH PC
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\radams\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\Robert\OneDrive\DATA\astra_toolbox\astra-1.8\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end

A00_prepareData
A01_getSourceData

% note, with paxscan according to varex:
% "will yield ~175µm real resolution with the detector’s 127 µm pixels"
% considering e.g. vibration probably reasonable to assume ~200 µm or worse in reality
% new sensor can have closer to 100 µm resolution, though, with Varex 1207

% first set of measurements:
% 8 mA beam current at 110 kV, spot size 1.2 mm
% 0.2 mm Cu filter, 0 mm Alu, "built-in" filter effectively 3.2 mm Alu @ 75 kV (??)
% source-detector 100 cm, obj-detector 5 cm
% about 1357 projections over 1 minute (fastest it could rotate @30 fps)
% ... using 1440 would be 0.25 degree steps\l
% with 0.1 mm resolution, ~30 mm FOV --> ~480 projections for full parallel sampling

pixelSizeFine = 0.01; % mm, pixel size for the finely defined phantom data
pixelSizeCoarse = 0.10; % mm, pixel size for image to be actually forward projected (must be multiple of pixelSizeFine!!)
sensorResolution = 0.2; % mm, assumed sensor resolution for tomography purposes
verticalResolution = 2; % assumed vertical binning of data to get better statistics, mm
nProjections = 1440; % rounded from above's 1357
%efficiency = 0.3; % wild guess for detector efficiency ... not used at the moment
beamCurrent = 8 ; % beam current in mA
totalMeasuringTime = 60; % total measuring time over all projections in seconds
projectionAngleListDeg = [360/nProjections:360/nProjections:360];

% 1 = with Film, 2 = empty, 3 = full
for objectCase = 1:2
% which Cu thickness to use for X-ray filtering...
%         case: 1   2   3   4   5   6   7   8
% Cu thickness: 0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0
for filterCase = 1:8
A04_prepareSourceData
% film case is film thickness * 0.1 mm, available points are 0.1:0.1:1.5 mm, or cases 1:15
for filmCase = 1:5
    
    % bins are from thisEnergy to thisEnergy+9 keV
    % max possible is 10:10:110 kV for 10-110 kV, i.e. 10:10:100
    % can also select sub-set of that
    energyBinList = 10:10:100;
    EbinCounter = 0; % for counting Ebin iterations
    for Ebin = energyBinList
        EbinCounter = EbinCounter+1;
        % this shifts the index by 9 due to the bin list starting at 10 keV
        Y_nps(EbinCounter) = sum(Y_sourceSpectrum(Ebin-9:Ebin*sourceBinSize)); % get #/pixel/projection total photons (assume 1 meter distance)
        disp(['Flat field photons per projection and pixel: ' num2str(round(Y_nps(EbinCounter))) ' in bin ' num2str(Ebin)])
        A03_preparePhantomData; % put xs data into phantom images in correct way
        % get diagonal distance in mm (via 2/sqrt(2) * side length),
        % and convert it to number of pixels but with a couple extras
        sensorPositions = 2*round(size(objectImageFilm,1)*pixelSizeCoarse/sensorResolution/sqrt(2)+2);
        % depending on which object case, use that object image
        if     objectCase==1; thisObjectImage = objectImageFilm;
        elseif objectCase==2; thisObjectImage = objectImageEmpty;
        elseif objectCase==3; thisObjectImage = objectImageFull;
        end
        % make astra volume geometry defining object space
        vol_geom = astra_create_vol_geom(size(thisObjectImage));
        % make astra projection geometry defining source-detector paths
        proj_geom = astra_create_proj_geom('parallel',...
            sensorResolution/pixelSizeCoarse, sensorPositions, projectionAngleListDeg*pi/180);
        proj_id = astra_create_projector('line',proj_geom, vol_geom);
        % put object image into astra form
        volume_id = astra_mex_data2d('create', '-vol', vol_geom, thisObjectImage);
        % do forward projection
        sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
        cfg = astra_struct('FP_CUDA');
        cfg.ProjectorId = proj_id;
        cfg.ProjectionDataId = sinogram_id;
        cfg.VolumeDataId = volume_id;
        fp_id = astra_mex_algorithm('create', cfg);
        astra_mex_algorithm('run', fp_id);
        disp('Forward projection complete.')
        % put non-noisy sinogram into a variable sinogramClean
        sinogramClean = astra_mex_data2d('get', sinogram_id);
        % get expected value of flat field counts at a given pixel
        EV(:,:,EbinCounter) = [Y_nps(EbinCounter)*exp(-sinogramClean)];
    end % end energy bin loop
    %
    for Ei = 1:EbinCounter % for each energy bin
        % add noise to the counts expected values
        noisyEV(:,:,Ei) = normrnd(EV(:,:,Ei),sqrt(EV(:,:,Ei)));
        % then weight that counts value by its energy
        % i.e. detector response is considered proportional to energy
        noisyEVweighted(:,:,Ei) = noisyEV(:,:,Ei) * energyBinList(Ei);
    end
    % sum the weighted noisy EV in energy from all bins
    sumNoisyEVweighted = sum(noisyEVweighted,3);
    % take the Y_nps and weight it as well by energy
    weightedFFcounts = dot(Y_nps,energyBinList);
    % calculate the noisy sinogram according to the noisy weighted values
    sinogramNoisy = -log(sumNoisyEVweighted/weightedFFcounts);
    % calculate noisy tomogram
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
    % garbage disposal
    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_projector('delete', proj_id);
    astra_mex_algorithm('delete', fbp_id);
    % save images
    if     objectCase==1; imgName = 'withFilm';
    elseif objectCase==2; imgName = 'empty';
    elseif objectCase==3; imgName = 'fill';
    end
    %
    if length(energyBinList)>1; spectrumType='polychromatic';
    elseif length(energyBinList)==1; spectrumType = ['monochromatic' num2str(energyBinList)];
    end
    %
    save(['./tomoData/tomo_filmCase_' num2str(filmCase) '_filterCase_' num2str(filterCase) '_' spectrumType '_' imgName '.mat'],'tomoNoisy')
    %
    
end % end film thickness case loop
end % end filter thickness case loop
end % end object case loop


