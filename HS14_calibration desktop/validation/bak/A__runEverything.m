

close all
clear
clc

% prepare astra toolbox stuff
% add mex and tools path for either personal laptop or ETH PC
thisDir = 'C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex';
A=exist(thisDir); if A==7, addpath(thisDir); end
thisDir = 'C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools';
A=exist(thisDir); if A==7, addpath(thisDir); end

addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign')
A00_prepareData


A01_getSourceData
%%
% spectrum & xs plot
globcol=get(groot,'DefaultAxesColorOrder'); %get colors
fh=figure(3249);clf
ax=gca;

%yyaxis left
ax.YScale='log';
hold on
grid on
name={'CHCl3','H2O','Alu'};
colind=[2,1,3]

for i =1:3 % attenuations
    p=semilogy(macroXSenergies,macroXSall(:,i),'Displayname',name{i},'color',globcol(colind(i),:));
    pub.BFLine(p,T.F)
end
ylh=ylabel('attenuation [1/cm]');
ylh.Color=[0 0 0];
pub.BFylab(ylh,T.F)

yyaxis right
ax=gca;
ax.YColor=[0 0 0]
plot(sourceDataAll(:,1),sourceDataAll(:,3)/10000000,'Displayname','x-ray spectrum','color',[0 0 0])
ylh=ylabel('normalized spectrum');

xlh=xlabel('Energy [keV]');
ylh.Color=[0 0 0];
lh=legend();
fname=sprintf('%sattenuation',T.Fig.saveto);

pub.BFfigure(fh,T.F)
pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)

savefig(fh,fname)
print(fh, '-dpdf',fname);
print(fh, '-dpng',fname);
clear fh ax name lh xlh ylh fname
%%
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
efficiency = 0.3; % wild guess for detector efficiency
beamCurrent = 8 ; % beam current in mA
totalMeasuringTime = 42; % total measuring time over all projections in seconds
projectionAngleListDeg = [360/nProjections:360/nProjections:360];

% 1 = with Film, 2 = empty, 3 = full
for objectCase = 1:2
    % which Cu thickness to use for X-ray filtering...
    %         case: 1   2   3   4   5   6   7   8
    % Cu thickness: 0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0
    for filterCase = 2
        A04_prepareSourceData
        % film case is film thickness * 0.1 mm, available points are 0.1:0.1:1.5 mm, or cases 1:15
        for filmCase = 3
            
            % bins are from thisEnergy to thisEnergy+9 keV
            % max possible is 10:10:110 kV for 10-110 kV, i.e. 10:10:100
            % can also select sub-set of that
            energyBinList = 10:10:100;
            EbinCounter = 0; % for counting Ebin iterations
            for Ebin = energyBinList
                EbinCounter = EbinCounter+1;
                % this shifts the index by 9 due to the bin list starting at 10 keV
                Y_nps(EbinCounter) = sum(Y_sourceSpectrum(Ebin-9:Ebin*sourceBinSize)); % get #/pixel/projection total photons (assume 1 meter distance)
                %disp(['Flat field photons per projection and pixel: ' num2str(round(Y_nps(EbinCounter))) ' in bin ' num2str(Ebin)])
                A03_preparePhantomData; % put xs data into phantom images in correct way
                % get diagonal distance in mm (via 2/sqrt(2) * side length),
                % and convert it to number of pixels but with a couple extras
                sensorPositions = 2*round(size(objectImageFilm,1)*pixelSizeCoarse/sensorResolution/sqrt(2)+2);
                % depending on which object case, use that object image
                if     objectCase==1; thisObjectImage = objectImageFilm;
                elseif objectCase==2; thisObjectImage = objectImageEmpty;
                elseif objectCase==3; thisObjectImage = objectImageFull;
                end
                %         % make astra volume geometry defining object space
                %         vol_geom = astra_create_vol_geom(size(thisObjectImage));
                %         % make astra projection geometry defining source-detector paths
                %         proj_geom = astra_create_proj_geom('parallel',...
                %             sensorResolution/pixelSizeCoarse, sensorPositions, projectionAngleListDeg*pi/180);
                %         proj_id = astra_create_projector('line',proj_geom, vol_geom);
                %         % put object image into astra form
                %         volume_id = astra_mex_data2d('create', '-vol', vol_geom, thisObjectImage);
                %         % do forward projection
                %         sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
                %         cfg = astra_struct('FP_CUDA');
                %         cfg.ProjectorId = proj_id;
                %         cfg.ProjectionDataId = sinogram_id;
                %         cfg.VolumeDataId = volume_id;
                %         fp_id = astra_mex_algorithm('create', cfg);
                %         astra_mex_algorithm('run', fp_id);
                %         %disp('Forward projection complete.')
                %         % put non-noisy sinogram into a variable sinogramClean
                %         sinogramClean = astra_mex_data2d('get', sinogram_id);
                n_ph = size(thisObjectImage,1);
                d_ph = pixelSizeCoarse*n_ph;
                p_d  = sensorResolution;
                n_d  = sensorPositions;
                sinogramClean=fanforward_new(thisObjectImage,...
                    n_ph,d_ph,p_d,n_d,projectionAngleListDeg*pi/180,950,50);
                
                % get expected value of flat field counts at a given pixel
                EV(:,:,EbinCounter) = [Y_nps(EbinCounter)*efficiency*exp(-sinogramClean)];
            end % end energy bin loop
            %
            for Ei = 1:EbinCounter % for each energy bin
                % add noise to the counts expected values
                noisyEV(:,:,Ei) = normrnd(EV(:,:,Ei),sqrt(EV(:,:,Ei)));
                noisyEVweighted(:,:,Ei) = noisyEV(:,:,Ei) * Y_nps(Ei); % weight it by its energy bin
            end
            % % % % % WARNING: MISSING AN ENERGY WEIGH TIN THE noisyEVweighted? I NEED TO CHECK THIS
            % sum the noisy EV in energy from all bins
            sumNoisyEVweighted = sum(noisyEVweighted,3);
            % get the flat field from all bins together weighted by energy
            weightedFFcounts = dot(Y_nps,energyBinList);
            % calculate the noisy sinogram according to the noisy weighted values
            sinogramNoisy = -log(sumNoisyEVweighted/weightedFFcounts);
            %             = -log(sum(noisyEV * Y_nps))
            % calculate noisy tomogram
            
            %     recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
            %     sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogramNoisy);
            %     cfg = astra_struct('FBP_CUDA');% myAlg
            %     cfg.ProjectorId = proj_id;
            %     cfg.ProjectionDataId = sinogram_id;
            %     cfg.ReconstructionDataId = recon_id;
            %     fbp_id = astra_mex_algorithm('create', cfg);
            %     astra_mex_algorithm('run', fbp_id);
            %     tomoNoisy = astra_mex_data2d('get', recon_id);
            %     disp('Filtered back projection complete.');
            %     % garbage disposal
            %     astra_mex_data2d('delete', sinogram_id, recon_id);
            %     astra_mex_projector('delete', proj_id);
            %     astra_mex_algorithm('delete', fbp_id);
            n_r=n_ph;
            p_r=d_ph;
            tomoNoisy=fan_new(sinogramNoisy,projectionAngleListDeg*pi/180,950,50,n_r,p_r,p_d);
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


function rec2=fan_new(sinogram,ang,src,det,n_r,p_r,p_d)
% sinogram= sinogram, size(sino) =[n_ang,n_pix]
% ang = angle list
% src = source - origin distance
% det = detector-origin distance
% n_r = width of recon volume in pixel
% p_r = recon volume pixel pitch in mm
%

d_r=p_r*n_r;
n_d=size(sinogram,2);
n_ang=length(ang);

vol_geom2 = astra_create_vol_geom(n_r, n_r,-d_r/2, d_r/2, -d_r/2, d_r/2);
proj_geom2 = astra_create_proj_geom('fanflat',p_d, n_d, ang, src, det);

%proj_id2 = astra_create_projector('strip_fanflat', proj_geom, vol_geom);
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom2, sinogram);
rec_id2 = astra_mex_data2d('create', '-vol', vol_geom2);

% create configuration
cfg2 = astra_struct('FBP_CUDA');
%cfg2.ProjectorId = proj_id2;
cfg2.ReconstructionDataId = rec_id2;
cfg2.ProjectionDataId = sinogram_id2;
%cfg.FilterType = 'Ram-Lak';
alg_id2 = astra_mex_algorithm('create', cfg2);
astra_mex_algorithm('run', alg_id2);

% get the reconstruction
%rec = astra_mex_data2d('get', rec_id2);
rec2 = astra_mex_data2d('get', rec_id2)*((p_r/p_d)^2);


rec2=rec2/p_d*n_ang/(pi/2); % <= Question 1) & 2)


astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');

end

function rec2=par_new(sinogram,ang,n_r,p_r,p_d)
% sinogram= sinogram, size(sino) =[n_ang,n_pix]
% ang = angle list
% src = source - origin distance
% det = detector-origin distance
% n_r = width of recon volume in pixel
% p_r = recon volume pixel pitch in mm
%

d_r=p_r*n_r;
n_d=size(sinogram,2);

vol_geom2 = astra_create_vol_geom(n_r, n_r,-d_r/2, d_r/2, -d_r/2, d_r/2);
proj_geom2 = astra_create_proj_geom('parallel',p_d, n_d, ang);

%proj_id2 = astra_create_projector('strip_fanflat', proj_geom, vol_geom);
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom2, sinogram);
rec_id2 = astra_mex_data2d('create', '-vol', vol_geom2);

% create configuration
cfg2 = astra_struct('FBP_CUDA');
%cfg2.ProjectorId = proj_id2;
cfg2.ReconstructionDataId = rec_id2;
cfg2.ProjectionDataId = sinogram_id2;
%cfg.FilterType = 'Ram-Lak';
alg_id2 = astra_mex_algorithm('create', cfg2);
astra_mex_algorithm('run', alg_id2);

% get the reconstruction
%rec = astra_mex_data2d('get', rec_id2);
rec2 = astra_mex_data2d('get', rec_id2)*((p_r/p_d)^2);
rec2 = rec2*p_d;  % fix 1 parallel

astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');

end

function sinogram=parforward_new(phan,n_ph,d_ph,p_d,n_d,ang)

% phan = phantom image
% n_ph = width of phantom in pixel
% d_ph = width of phantom in mm
% p_d = detector pixel pitch
% n_d = detector pixel amount
% and = angle list in radian

proj_geom = astra_create_proj_geom('parallel', p_d ,n_d, ang);
vol_geom = astra_create_vol_geom(n_ph,n_ph,-d_ph/2,d_ph/2,-d_ph/2,d_ph/2);

% store volume
volume_id = astra_mex_data2d('create', '-vol', vol_geom, phan);

% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id);

sinogram=sinogram/p_d;  % fix 1 parallel

% garbage disposal
astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');
end

function sinogram=fanforward_new(phan,n_ph,d_ph,p_d,n_d,ang,src,det)

% phan = phantom image
% n_ph = width of phantom in pixel
% d_ph = width of phantom in mm
% p_d = detector pixel pitch
% n_d = detector pixel amount
% and = angle list in radian

proj_geom = astra_create_proj_geom('fanflat',p_d, n_d, ang, src, det);
vol_geom = astra_create_vol_geom(n_ph,n_ph,-d_ph/2,d_ph/2,-d_ph/2,d_ph/2);

% store volume
volume_id = astra_mex_data2d('create', '-vol', vol_geom, phan);

% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id);

% garbage disposal
astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');
end

