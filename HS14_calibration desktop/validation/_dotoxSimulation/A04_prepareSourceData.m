

% calculated values
exposurePerProjection = beamCurrent*totalMeasuringTime/nProjections;
disp(['Exposure per projection = ' num2str(exposurePerProjection)]);
pixelSizeRatio = pixelSizeCoarse/pixelSizeFine; % ratio of coarse to fine
if floor(pixelSizeRatio)~=(pixelSizeRatio),disp('Error! PixelSizeRatio is non-integer.');return;end
% source energy steps of 1 keV from 10 to 110
E_sourceSpectrum = sourceDataAll(:,1);
% source output photons/keV/cm^2/mAs @ 1 meter
Y_sourceSpectrum = sourceDataAll(:,filterCase+1);
sourceBinSize = E_sourceSpectrum(2)-E_sourceSpectrum(1); % energy bin of source spectrumd data, keV
% the scripts assume 1 keV bins, if not then everything needs to be checked:
if sourceBinSize~=(1),disp('Error! Source bin size not 1 keV as expected. Check everything.');return;end
if pixelSizeFine~=(0.01),disp('Error! pixelSizeFine not = 0.01 mm as expected. Check everything.');return;end
Y_sourceSpectrum = Y_sourceSpectrum * verticalResolution/10 * sensorResolution/10 * exposurePerProjection ; % convert to #/keV/pixel/projection @ 1 meter