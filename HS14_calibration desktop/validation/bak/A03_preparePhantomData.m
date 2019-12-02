


    % vectors already defined at 10:1:110 kV [1/cm]:
    % macroXSenergies
    % macroXSaluminum
    % macroXSchloroform
    % macroXSwater
    thisXSaluminum = mean(macroXSaluminum(Ebin-9:Ebin*sourceBinSize));
    thisXSchloroform = mean(macroXSchloroform(Ebin-9:Ebin*sourceBinSize));
    thisXSwater = mean(macroXSwater(Ebin-9:Ebin*sourceBinSize));
    
    % phantom raw data:
    % pixel resolution 0.01 mm
    % 0 = aluminum, 1 = outside channel, 2 = H2O, 3 = void in channel, 4 = film
    % prepare image and apply attenuation values:
    objectImageFilm = phantomDataAll(:,:,filmCase);
    % define the three regions which don't change
    objectImageFilm(objectImageFilm == 0 ) = thisXSaluminum;
    objectImageFilm(objectImageFilm == 1 ) = 0;
    objectImageFilm(objectImageFilm == 2 ) = thisXSwater;
    % also make empty and full versions
    objectImageEmpty = objectImageFilm;
    objectImageFull = objectImageFilm;
    % real image with finite film
    objectImageFilm(objectImageFilm == 3 ) = thisXSchloroform*4.475/1411; %cl3 vap;
    objectImageFilm(objectImageFilm == 4 ) = thisXSchloroform;
    % empty
    objectImageEmpty(objectImageEmpty == 3 ) = thisXSchloroform*4.475/1411; %cl3 vap;;
    objectImageEmpty(objectImageEmpty == 4 ) = thisXSchloroform*4.475/1411; %cl3 vap;;
    % full
    objectImageFull(objectImageFull == 3 ) = thisXSchloroform;
    objectImageFull(objectImageFull == 4 ) = thisXSchloroform;
    
    objectImageFilm = objectImageFilm * pixelSizeCoarse/10; % convert 1/cm to 1/pixel
    objectImageEmpty = objectImageEmpty * pixelSizeCoarse/10; % convert 1/cm to 1/pixel
    objectImageFull = objectImageFull * pixelSizeCoarse/10; % convert 1/cm to 1/pixel
    
    objectImageFilm = imresize(objectImageFilm,1/pixelSizeRatio);
    objectImageEmpty = imresize(objectImageEmpty,1/pixelSizeRatio);
    objectImageFull = imresize(objectImageFull,1/pixelSizeRatio);