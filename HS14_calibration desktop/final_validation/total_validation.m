% total validation scrip according to roberts ideas

close all
clearvars -except T F draw d mtomo ratio
clc

addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\apriori')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final_validation')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\validation\_dotoxSimulation')

% load all required data


% load real data
if 0
    load('T'); % get recon data
    load('F'); % get lft data
    cas=1;rep=1;T.cas=cas;T.rep=rep;
    
    %raw data
    dRaw=f.ReadXrayTomo(T);
    fname=sprintf('%s%02d_%02d.mat','1_corr_',T.cas,T.rep);
    fpath=sprintf('%s%s',T.d.DataPath,fname);
    
    % processed data
    d=f.BoLoad(fpath,T);

    %miachas tomos
    %mtomo=f.ReadMichasTomos(F,T.Rec.recsize,210);
    load('mtomo')
    load('ratio')
end


%%
% parameters

geo.ph_n_coarse = 0.127; % mm, pixel size for image to be actually forward projected (must be multiple of pixelSizeFine!!)
geo.phanResizeRatio=11;
geo.ph_n_fine = geo.ph_n_coarse/geo.phanResizeRatio; % mm, pixel size for the finely defined phantom data
geo.det_p = 0.127; % mm, assumed sensor resolution for tomography purposes
geo.det_vert = 2; % assumed vertical binning of data to get better statistics, mm
geo.nProjections = 1357; % rounded from above's 1357
geo.efficiency = 0.3; % wild guess for detector efficiency
geo.beamCurrent = 8 ; % beam current in mA
geo.totalMeasuringTime = 45.120; % total measuring time over all projections in seconds
geo.ang=deg2rad(T.Rec.angles);
geo.n_d=269; % number of detector pixels used
geo.src=950; % source origin distance in mm
geo.det=50;% source detector distance in mm



% energy & att coefficients
[energy,att]=val.get_xs(); %in keV and 1/cm

%phantoms
lftStep=0.05; % film thickness step size & minimal thickness
lftMax=0.5; % maximal film thickness
lfts=0:lftStep:lftMax; % lft set
width=geo.n_d*geo.phanResizeRatio; %phantom size in pixels
if 0 % put 0 to save time
    [phan,phanCoord]=val.makePhantomSeries(lfts,width,geo.ph_n_fine);
    save('phan','phan')
    save('phanCoord','phanCoord')
else
    load('phan')
    load('phanCoord')
end
val.plotPhans(phan,phanCoord,T);

%spectrum
[spec,cufilts]=val.get_spectrum(); % in N[keV cm^2 mAs]^-1 @ 1 meter

% plot attenuation & spectrum
val.plotxs(energy,att,spec,T);
val.plotspec(energy,spec,cufilts,T)


%%

% get the pixel flux
pixelFlux=val.getFluxSet(geo,spec,energy,cufilts); % in [N/keV/pixel] @ 1m 

% define the energy binning
eBinning=[10];% [1,5,10,20]; %number of energy values in one bin, array
tic
clear binCenter binRange specbin fluxbin rec sinoClean
for i=1:length(eBinning)
    energyBinSize=eBinning(i); % number of energy value, single integer
    neBins=100/energyBinSize; % number of bins
    energyBinStart = 11:energyBinSize:110; %energy(1):engergy(end)
    energyBinInd=energyBinSize:energyBinSize:100;
    specbin=zeros(length(energyBinStart),length(cufilts));
    attbin=zeros(length(energyBinStart),size(att,2));
    
    % make energy binned properties
    for en=1:length(energyBinStart)
        %disp(en)
        % set sum range and sum over it
        binrange=(energyBinInd(en)-(energyBinSize-1)):energyBinInd(en);
        binCenter(en)=mean(energy(binrange));
        
        % binned XS
        for mat=1:size(att,2)
            attbin(en,mat)=mean(att(binrange,mat)); % 1/cm
        end
        
        % binned Spectrum
        for filter =1:length(cufilts)
            specbin(en,filter)=sum(spec(binrange,filter));% N[keV cm^2 mAs]^-1 @ 1 meter
            fluxbin(en,filter)=sum(pixelFlux(binrange,filter)); % Roberts Y_nps in [N/keV/pixel] @ 1m 

        end
    end
    %plot bins
    val.plotBinningAtt([1,1,1],energy,att,attbin,energyBinStart,energyBinSize,T)
    val.plotBinningSpec([1,1,1],energy,spec,specbin,attbin,energyBinStart,energyBinSize,cufilts,T)
    

    % apply phantom attenuations
    attphan=val.ImbuePhantom(phan,attbin/10,geo); % factor 10 from 1/cm to 1/mm
    
    % forward projection
    sinoClean=zeros(geo.nProjections,geo.n_d,length(lfts),neBins); % in [-],  1/cm*cm
    for lft =1:length(lfts)
        %disp(lft)
        for en=1:length(energyBinStart)
            sinoClean(:,:,lft,en)=val.fanforward_new(attphan(:,:,lft,en),...
                size(attphan,1)*geo.ph_n_coarse,geo.det_p,geo.n_d,...
                geo.ang,geo.src,geo.det);
        end
    end
    sinoClean=convn(sinoClean,mt.ElliKernelGen(1,3,1,1),'same'); %add blurring, arbitrarily chosen
    %
    for filter =[2];%1:length(cufilts)
        
        % make a temporary flux array
        fa=permute(repmat(fluxbin(:,filter),...
            [1,geo.nProjections,geo.n_d,length(lfts),1]),[2 3 4 1]);
        
        % also make a temporary weight array, reflecting the
        % assumption that the detector pixel count for a photon is
        % proportional to its energy
        wa=permute(repmat(binCenter'/mean(binCenter),...
            [1,geo.nProjections,geo.n_d,length(lfts),1]),[2 3 4 1]);
        % expected pixel values as efficiency*e(-sino)*flux per energy bin
        % for each filter in [N/keV/pixel] @ 1m 
        EVfilt=exp(-sinoClean).*fa;
        
        % add normal distributed noise
        EVfiltNoised=normrnd(EVfilt,sqrt(EVfilt));
        
        % add the energy weighting
        EVfiltNoisedWeighted=EVfiltNoised.*wa;

        % make the flat field
        % skipped
        sino=-log(sum(EVfiltNoisedWeighted,4)./sum(fa.*wa,4));
        %             for filter =1:length(cufilts)
        %                 sinoNWeighted(:,:,lft,en)=sinoNoised(:,:,lft,en)*fluxbin(:,:,en,filter)
        %             end
        
        % it's recon time!
        for lft=1:length(lfts)
            rec(:,:,lft,filter)=val.fan_new(sino(:,:,lft),geo.ang,...
                geo.src,geo.det,geo.n_d,geo.ph_n_coarse,geo.ph_n_coarse);
        end
        fid=10000+filter+10*i;
        val.HistCompare(fid,sino(:,:,1),...
            squeeze(d(:,100,:))',...
            cufilts(filter),neBins,T)
        
        phanweights=permute(repmat(fluxbin(:,filter)/mean(fluxbin(:,filter)),...
            [1,size(attphan,1),size(attphan,2),length(lfts),1]),[2 3 4 1]);
        phanweighted=mean(attphan.*phanweights,4);
        %phanrange=16:284;
        
        % diff recon plots
        val.ShowRecons(rec(:,:,:,filter),lfts)
   
    end
   val.ReconHistCompare(phanweighted(:,:,1),...
            rec(:,:,1,filter),mtomo(:,:,1,8),cufilts(filter))
        
    toc
end

%% do the lft analysis




phanprof=zeros(F.pth.r_path+1,length(cufilts),length(lfts),F.geo.n_pins,F.pth.n_angles);
phandifprof=zeros(F.pth.r_path+1,length(cufilts),length(lfts)-1,F.geo.n_pins,F.pth.n_angles);


for filt=2 %length(cufilts)
    for lft=1:length(lfts)
        im=rec(:,:,lft,filt);
        for pin=1:F.geo.n_pins
            xstart=F.cen.cen(pin,1);
            ystart=F.cen.cen(pin,2);
            for ang=1:F.pth.n_angles
                xstop=F.Pa.ProfEndPnt(pin,ang,1);
                ystop=F.Pa.ProfEndPnt(pin,ang,2);
                phanprof(:,filt,lft,pin,ang)=...
                    improfile(im',[xstart,xstop],[ystart,ystop],F.pth.r_path+1);
                
            end
        end
        if lft >1
            phandifprof(:,filt,lft-1,:,:,:)=phanprof(:,filt,lft,:,:,:)-phanprof(:,filt,1,:,:,:);
        end
    end
end

%% check the profs
angs=1:10:F.pth.n_angles;
nang=length(angs);
filt=2;
lft=1;
clear x y plotprof
for pin=1:4
    for ang2=1:nang;
        ang=angs(ang2);
        x(:,pin,ang2)=linspace(F.cen.cen(pin,1),F.Pa.ProfEndPnt(pin,ang,1),size(F.difProf,5));
        y(:,pin,ang2)=linspace(F.cen.cen(pin,2),F.Pa.ProfEndPnt(pin,ang,2),size(F.difProf,5));
        plotprof(:,pin,ang2)=phanprof(:,filt,lft,pin,ang);
    end
end
[fh,ax]=pub.checkprofile(im,plotprof,x,y);
%% 
figure(23462);clf
pin=1;
col=hsv(nang);
for lft=2:length(lfts)
    subplot(1,length(lfts)-1,lft-1)
    hold on
    for ang=1:nang
        plot((1:size(F.difProf,5))*F.recres,squeeze(phandifprof(:,filt,lft-1,pin,ang)),...
            'color',col(ang,:),'displayname',sprintf('%d°',angs(ang)));
    end
    ylim(0.06*[-0.1,1])
    xlabel('path length [mm]')
    title(sprintf('sim %3.2fmm lft',lfts(lft)))
    grid on
    
    xlim([0 7])
    
    if lft==2
        legend();
        ylabel('\mu [1/mm]')
    end
end
%% plot the line integrals
% first sum up the profiles
phanmusraw=squeeze(sum(phanprof(F.pth.lftSumRange,:,:,:,:),1));

%then make the difference
clear phanmus
for i =1:(length(lfts)-1)
    phanmus(:,i,:,:)=phanmusraw(:,i+1,:,:)-phanmusraw(:,1,:,:);
end

% now average over angles
phanmusangmean=mean(phanmus(:,:,:,F.pth.angrang),4);
%% now plot
figure(983);clf;
pin=1;
plot(lfts,[0,squeeze(phanmusangmean(filt,:,pin))],'-x')
xlabel('lft [mm]')
ylabel('\Sigma \mu(s)\cdots  [-]')
grid on
title('[sim] path line integral vs. film thickness')

figure(984);clf;
pin=1;
plot(lfts,[0,squeeze(phanmusangmean(filt,:,pin))]./lfts,'-x')
xlabel('lft [mm]')
ylabel('\mu_{eff} [1/mm]')
grid on
title('[sim] resulting effective attenuation')


%% phantom size study
[tph,tphC]=val.makePhantomSeries(0,width,geo.ph_n_fine);
[ph2,Co22]=val.make_phantom(0,11*269,0.127/11);
co3=imresize(tphC(:,:,1),1/11,'box');
co3(1,1)+0.127*269/2

%%
res1=0.127;         % resolution 1 or pixel width
npix1=269;       % number of pixel withd of image
numxmax1=(npix1-1)/2; % maximum pixel number in each direction
numx1=linspace(-numxmax1,numxmax1,npix1); % number of each pixel
coxmax1=res1*(npix1-1)/2; % coordinate of the furthest pixel
cox1=linspace(-coxmax1,coxmax1,npix1); % coordinate list

upsample=11;    % upsampling factor

res2=res1/upsample;     % upsampled resolution
npix2=npix1*upsample;   % number of pixels when upsampled
numxmax2=(npix2-1)/2;   % maximum pixel number in each direction
numx2=linspace(-numxmax2,numxmax2,npix2); % number of each pixel
coxmax2=res2*(npix2-1)/2; % coordinate of the furthest pixel
cox2=linspace(-coxmax2,coxmax2,npix2); % coordinate list usampled

[xx,yy]=ndgrid(cox2,cox2);
xxx=imresize(xx,[npix1,npix1],'box','Method','bilinear','Antialiasing',false);


%%
xmax=10;
weite=2*10+1;
sample=11;

xmax2=xmax+10*0.127/11
x=linspace(-xmax2,xmax2,269*11)/2;
x=linspace(-(10+5/11),(10+5/11),121)/2;
[xx,yy]=ndgrid(x,x);
xxx=imresize(xx,[11,11],'box','Method','bilinear','Antialiasing',false);
xxx(1,1)
xmax/2
mean(mean(xx(1:11,1:11)))
%%
% 1) make upsampled phantom
phan3=phan(:,:,1);phan3(phan3==1)=3;phan3(phan3==2)=3;phan3(phan3==0)=1;phan3(phan3==3)=0;
phan4=imresize(phan3,1/11);
phan42=imresize(phan3,1/11,'box');
figure(23451);clf
imshow(phan4-phan42,[])
figure(26743);clf;
ax(1)=subplot(1,3,1);
imshow(phan4,[]);
ax(2)=subplot(1,3,2);
imshow(phan42,[]);
ax(3)=subplot(1,3,3);
imshow(phan42-phan4,[]);
linkaxes(ax)



% 2) make regular phantom


% 3 compare

    
%% make phantoms for micha
% 1) analytic

Ph=f42_RotPhan(f42_ChannelPhanTransform,pi/4);
li=Ph.CL==0;
ci=Ph.CL==1;
MichaPhanLine = table(...
    (1:sum(li))',... % counter
    Ph.reg(li,1),...
    Ph.reg(li,1),...
    Ph.Lx1(li,1),...
    Ph.Lx1(li,2),...
    Ph.Lx2(li,1),...
    Ph.Lx2(li,2),...
    'VariableNames',{...
    'Nr',...
    'Region_swich_1',...
    'Region_swich_2',...
    'Line_start_XY_1',...
    'Line_start_XY_2',...
    'Line_end_XY_1',...
    'Line_end_XY_2'})

MichaPhanCirc=table(...
    (1:sum(ci))',...
    Ph.c(ci,1),...
    Ph.c(ci,2),...
    Ph.r(ci)',...
    Ph.p0(ci)',...
    Ph.dp(ci)',...
    Ph.reg(ci,1),...
    Ph.reg(ci,1),...
    'VariableNames',{...
    'Nr',...
    'center_XY_1',...
    'center_XY_2',...
    'Radius',...
    'Start_Angle',...
    'Angular_Distance',...
    'Region_swich_1',...
    'Region_swich_2'})
writetable(MichaPhanCirc,sprintf('%sPhantom_Circles.txt',T.Rec.michamainpath)...
    ,'Delimiter','space')
writetable(MichaPhanLine,sprintf('%sPhantom_Lines.txt',T.Rec.michamainpath),...
    'Delimiter','space')


% coordinate list
res1=0.127;         % resolution 1 or pixel width
npix1=269;       % number of pixel withd of image
numxmax1=(npix1-1)/2; % maximum pixel number in each direction
numx1=linspace(-numxmax1,numxmax1,npix1); % number of each pixel
coxmax1=res1*(npix1-1)/2; % coordinate of the furthest pixel
cox1=linspace(-coxmax1,coxmax1,npix1); % coordinate list
coordtable=table((1:269)',cox1','VariableNames',...
    {'pixel_number','pixel_center_coordinate'});

writetable(coordtable,sprintf('%sPixel_Coordinates.txt',T.Rec.michamainpath),...
    'Delimiter','space')












