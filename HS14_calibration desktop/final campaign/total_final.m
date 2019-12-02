%%

% this is the script that evaluates the data from 
% the final measurement campaign in Descember 2017

% data is being held on disk, parameters are saved in struct T (Total)
% T total struct, frequently used parameters go here
% T.d data (raw) related properties
% T.dose contains dose corection related stuff
clc
input('really runn the total script? abort now^^');

%% initialization
clc             % clear commad line output
clear global    % clear globals
clear           % clear workspace
close all       % close figures
format compact  % compact command line output
T.sys.t_tot=tic;
T.sys.StoreFormat='uint16'; % how arrays shall be sotred
tic;
fprintf('running the tomo processing script...\n')

% .d stands for data, raw data, and related properties
T.d.DataPath='E:\20171213 final campaign\'; % dir of raw data

%T.d.DataPath='C:\data\final campaign\';
T.proofs=1; % flag to activate (or deactivate) proloadoving routines & user interactivities
T.d.header=2048; % the file header to be skipped at read in
T.d.imsize=[640 1024]; % the frame size in pixels
T.d.GreyBytes=2; % the color depth in bytes, usually 2 for VIVA .seq videos
T.d.ncas=9; %number of cases/ Betriebspunkte
T.d.nrep=10; % number of repetitions for each case

% generate list of raw data files
T.d.fileflag=('*tom_*.seq'); % files with this flag will enter the file list
[T.d.List,T.d.nframes,T.d.reps]=f.GenerateFileList(T,10);    % finally read in the raw files list

% pre-readin & generate sino, rad & yshift-plane
T.Raw.sinoyplanes=[50];% yplanes for raw (angular gauge), bottom and up sinogram

% crop parameters
T.Raw.CropThrash=5000; % guess Value for the "find the rotation start & end" algo
T.Raw.CropMargin=10; % number of pixels to be kept before and after the crop
T.Raw.RotCutLimit=1460; % within how many pixels shall we look for the end 
                        % of rotation? =pixel length of rotation + ~20

% dose correction & beam fix
% specify the regions that the dose calculation is based on:
T.dose.windows=[5 5 455;... %x start
    200 200 635;...         %x end
    110 800 110;...         %y start
    325 1018 1018];         %y end

% matching settings
T.Match.SinoWindows=[5,455;... % x start
        195,635;            % x end
        20,20;...         % frame start
        1443,1443];           % frame end
T.Match.GuessLine=150; %x-value for the xcov guess 

T.Match.YWindows=[136;... % y start
        348;            % y end
        30;...         % frame start
        1460];           % frame end
T.Match.yshiftplane=322;% x-value, at which the y shift plane lies
T.Match.YGuess=2;   % guess value for that shift. shift will be searched in
T.Match.YRange=4;      % a Yshift search range (+-)
T.Match.YSteps=21;      % number of steps that cover +-Yrange
T.Match.PinMaskFlag=1; % should the experimental pin mask be used?
T.Match.MaskThresh=0.6; % max. grey value for the pin mask
T.Match.DilBox=5;       % dilation box: pin mask will be enlarged by it

% the quest for 360°
T.q360.startframe=80; % this is will be the first frame of the 360° in the 
                      % cropped sinogram. ideally where the angular gauges
                      % are both +-45° relative to the beam (works with
                      % less though)
T.q360.NFGuess=1357;  % guess of the number of frames for 360°
T.q360.NFGrange=50;   % range to look within for the 360°
T.F=f.FigProperties(T);   % centralized figure properties                  



% centering & reconstruction
T.Cen.start=193; % rough cut to center, in terms of x-pixel
T.Cen.stop=461;
T.Cen.range=T.Cen.start:T.Cen.stop; % in here lies the channel
T.Rec.recsize=length(T.Cen.range); % reconstruction edge size
T.Rec.angles=linspace(0,360,T.q360.NFGuess+1); %make the angles list
T.Rec.angles(end)=[];
T.Cen.CSFR{2}=345:1010; % CentShiftFitRange: y-planes, based on which the 
T.Cen.CSFR{1}=400:910;  % center shift will be performed
T.Rec.SinoRotShift=200; % sino is circshifted by this to have the 4 quadrants
                        % aligned
T.Rec.nRcas=4; % number of ReconCases
T.Rec.R2Pcas=[1 6 7 8]; % ReconCas to Preprocesing cas
T.Rec.Rcas=[1 2 3 4 5 2 3 4 9 ]; % translates the cas 1678 to 1234;
R2P=T.Rec.R2Pcas;                        
T.Cen.quadrot=70; % rotate to have the channel quadrants aligned
T.Cen.xmin=-2;
T.Cen.xmax=2;
T.Cen.MaskThresh=0.003;
T.Cen.FiboThresh=0.007;
T.Cen.MaskPlane=100; % 
T.Cen.niter=20;

T.Rec.detPitch=0.127; % detector pixel pitch
T.Rec.src=950;          % source axis distance in mm
T.Rec.det=50;           % detector axis distance in mm
T.Rec.michamainpath='T:\cbolesch\Sino_fuer_micha\';


%preallocate, CAUTION deletes existing values
T=f.PreProcessingPreAllocation(T);
%%

% parforhacks
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

% Preprocessing
tic
for cas=1:T.d.ncas %iterate cases
    if any(T.d.nframes(cas,:)~=0) % skip if no files available
    T.cas=cas; % copy iterater to T struct
    T.Rec.michapath{cas}=sprintf('%s%d\\',T.Rec.michamainpath,cas);
    mkdir(T.Rec.michapath{cas});
    for rep=1:T.d.nrep % iterate repetitions
        d=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
        if T.d.nframes(cas,rep)~=0 % skip if no measurement available

            T.rep=rep; % copy iterater to T struct
            fprintf('%s preprocessing...\n',f.CommLineStart(cas,rep)); %debug
            if 0 % skip some data processing steps if we did them already
                % 1 read in
                d=f.ReadXrayTomo(T);
                T.Raw.rad(:,:,cas,rep)=squeeze(d(:,:,1));
                
                % 2 find rot startt & stop and crop
                T.Raw.CropRange(cas,rep,:,:)=f.FindCrop(T,...
                    squeeze(d(:,T.Raw.sinoyplanes(1),:)));
                
                % make 360° check plot
                T.q360.nFrames(cas,rep)=...
                    f.frames360(squeeze(d(:,T.Raw.sinoyplanes(1),:)),T);
                
                % plot the graph
                f.CheckSino1(T,d,10);
                
                % crop & copy
                [d,T.Raw.nframes(cas,rep)]=f.SinoCrop(d,T);
                
                % 3 dose correction, some plots
                [d,T.dose.doses(cas,rep,:),dosemin]=f.DoseCorrect(d,T);
                f.DoseMaskCheck(d(:,:,1),T,40);
                f.DoseCorrectCheck1(T,dosemin,50);
                f.DoseCorrectCheck2(d,T,70);
                clear dosemin
                
                % 4 beam correct
                
                [d,T.fit.fits{cas,rep},T.fit.fitdose{cas,rep},im_old]=...
                    f.BeamCorrect(d,T);
                f.BeamCorrectCheck(d,T,im_old);
                f.DoseCorrectCheck2(d,T,71);
                clear im_old
            else
                fname=sprintf('%s%02d_%02d.mat','1_corr_',T.cas,T.rep);
                fpath=sprintf('%s%s',T.d.DataPath,fname);
                % temporarily change it to 
                fpath=sprintf('%s%s%s',T.d.DataPath,'back 4rec\',fname);
                d=f.BoLoad(fpath,T);
            end
            
            % crop z
            d=d(:,:,T.q360.startframe:(...
                T.q360.startframe+T.q360.nFrames(1,1)-1));
            
            % image correction: pixels & filter
            corr=f.ImCorrFilt(mean(d,3),T);
            % apply image correction
            d=f.ApplyImCorr(d,corr,T);
            
            % crop x & y
            d=d(T.Cen.range,:,:);
            
            % quadrant-rotate, fliplr (parity), log
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
            
            v=nan(niter,BS(2));
            xlist=nan(niter,BS(2));
            eshift=nan(1,BS(2));
            recon=zeros(recsize,recsize,BS(2),'single');

                        
            % files for micha

                            % save intermediate result
            T.fnames.corr{cas,rep}=f.BoSave3(d,'1_corr_',T);
            
            if 1 % temporarily use the previous shift values
            eshift=T.Cen.eshift(:,cas,rep);
            xlist=T.Cen.xlist(:,:,cas,rep);
            v=T.Cen.v(:,:,cas,rep);
            end
            
            fprintf('planes: ')
            for plane=10:1017% plane=[50,512,1000]% 
                
                f.f_BoCount(plane,100,10,5)
                try
                    % check out sinogram
                    sino=squeeze(d(:,plane,:));
                    if 0
                        % center exact
                        [eshift(plane),xlist(:,plane),...
                            v(:,plane)]=f.fibo(sino,0,...
                            xmin,xmax,niter,0.007,1,ang,recsize,detPitch,src,det,MMask);
                        eshift(plane)=round(eshift(plane),2);
                        
                    end
                    
                    % reconstruct
                    recon(:,:,plane)=val.fan_new(...
                        f.fraccircshift(sino,-eshift(plane))',...
                        ang,src,det,recsize,T.Rec.detPitch,T.Rec.detPitch);
                        
                    
%                     recon(:,:,plane)a.FBPexplFan(...
%                         f.fraccircshift(sino,-eshift(plane))',...
%                         recsize,ang,detPitch,src,det);
                    

%                     f.CheckRecon(recon(:,:,plane)...
%                         ,v(:,plane),xlist(:,plane)...
%                         ,cas,rep,plane,xmin,xmax,1)
                    
                    %save to disk
                    fout=sprintf('%sMicha_c%d_r%d_p%04d.csv',T.Rec.michapath{cas},cas,rep,plane);
                    dlmwrite(fout,recon(:,:,plane),'delimiter',' ','precision',7);
                    
                    
                catch e %e is an MException struct
                    fprintf(2,'%s\n',e.identifier);
                    fprintf(2,'%s\n',e.message);
                    fprintf('%4d didnt reconstruct\n',plane)
                end
            end
            fprintf('\n')
            T.Cen.eshift(:,cas,rep)=eshift;
            T.Cen.xlist(:,:,cas,rep)=xlist;
            T.Cen.v(:,:,cas,rep)=v;
            T.fnames.rec{cas,rep}=f.BoSave3(recon,'4_rec_',T);

            
            
            
        end   % isempty(T.fnames.add{cas}) % skip if no files available
    end
    end
end
save('T','T')
%% post process the micha-matches
% F the film struct
clear 
clc
close all
load('T')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\apriori')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final_validation')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\validation\_dotoxSimulation')
F.pp.loadpath='T:\cbolesch\res'; % the path for loading the recons
F.pp.casstr=8;
F.pp.plastr=11:14;
F.pp.description='stuff related to reading in the micha files';
[F.pp.flist,F.pp.cases,F.pp.planes,F.raw]=f.GenRecFileList(F,200);
F.planes=10:16:1002; % manual planes fix

F.raw.maxp=max(F.raw.nplane);
mtomo=f.ReadMichasTomos(F,T.Rec.recsize,210);

%
F.pth.description='parameters involved in the profile paths';
F.pth.r_path=60; % path length in pixel +manual fix

F.size=size(mtomo);
F.Ncases=F.size(3);%F.size(1);
F.recsize=F.size(1);
F.h=F.size(4); % stack height
F.recres=0.127; % resolution of the reconstructed image

F.geo.description='geometrical parameters';
F.geo.r_rod=5.14;       % rod outer radius
F.geo.dh=0.8473;         % hydraulic diameter in cm
F.geo.n_pins = 4; % number of fuel pins
F.geo.realangle=174.83; % manual check in inventor

F.pth.lftSumRange=38:46; % empirical best fit
F.pth.d_angle=180; % angle range
F.pth.angrang=4:177; % indices of valid angles
F.pth.n_angles=F.pth.d_angle+1; % number of angles

F.par.flow=[0 3.45,2.87,2.3];
F.par.flowkgs=[0 5.00,4.05,3.25];
F.par.RobXS=0.059; % robert's 'magic' blackbox cross section, in 1/mm
F.XS=F.par.RobXS;
%F.att=xs; %function xs for attenuation coefficients

F.pth.a0=deg2rad([-90,-180,90,0]); % first angle for each circle
F.pth.Profs=zeros(F.Ncases,length(F.planes),F.geo.n_pins,F.pth.n_angles,F.pth.r_path+1);


% y coordinates in hydraulic diameters:
% plane(22) is the first plane downstream of the spacer where michas algo
% converges
F.par.description='random parameters of all sorts';
F.par.spacerplane=292; % the plane that the spacer ends; zero-reference

%create the imaging coordinate system, plot parameters etc
F.p.description='parameters for plotting, i.e. scales of hydraulic diameter, mm, planes, coordinates,...';
F.p.dh=(F.planes-F.planes(22))*T.Rec.detPitch/(F.geo.dh*10); 
F.p.xcoord=lft.coords(F); % tomo pixel coords
F.p.RA=imref2d([F.recsize,F.recsize],F.p.xcoord([1,end]),F.p.xcoord([1,end]));
F.p.mPlanes=1:F.h; % number of michas planes
F.p.mPlanesOrig=F.p.mPlanes*16+10; % each first planes over which micha averaged
F.p.mPlanesmm=F.p.mPlanesOrig-F.par.spacerplane*T.Rec.detPitch; % mm relative to spacer edge
F.p.dhticklabels=[-3:3:9]; % d_h positions to put labels on in plots
F.p.dhticks=(lft.dh2plane(F.p.dhticklabels)-10)/16; % 
%
fh=figure(312);clf
imshow( mtomo(:,:,1,50)',[])
set(gca,'YDir','normal')
hold on
colorbar
axis on
%
F.cen.prot=45; % deg
F.cen.pzoom=1/F.recres; % factor
F.cen.px=ceil(F.recsize/2);
F.cen.py=ceil(F.recsize/2);
F.cen.pcenid=[61 13 29 45]; % counter-clockwise, starting east


F.cen.P=f42_ChannelPhanTransform();
F.cen.P2=f42_RotPhan(F.cen.P,deg2rad(F.cen.prot));
F.cen.P2=f42_PhanZoom(F.cen.P2,F.cen.pzoom);
F.cen.P2=f42_PhanMove(F.cen.P2,F.cen.px,F.cen.py);
f42_PlotPhan(F.cen.P2,gca,[1 0 0],[1 4]);

F.cen.cen=F.cen.P2.c(F.cen.pcenid,:);
%
for i =F.cen.pcenid%1:length(P2.c)
    plot(F.cen.P2.c(i,1),F.cen.P2.c(i,2),'ob');
    %text(F.cen.P2.c(i,1),F.cen.P2.c(i,2),num2str(F.cen.P2.r(i)),'FontSize',20);
end
%
% find paths
[F.Pa.ProfEndPnt,F.Pa.ProfAngles]=lft.FindProfilePathsM(F);

% check path enpoints
for i=1:4
    plot(squeeze(F.Pa.ProfEndPnt(i,:,1)),squeeze(F.Pa.ProfEndPnt(i,:,2)),'.m')
end

save('F','F')


% obtain the profiles
for i = 1:F.raw.ncas
    F.pth.Profs(i,:,:,:,:)=lft.FindProfilesM(F,squeeze(mtomo(:,:,i,:)));
end
% make differential profs
for i=1:(F.raw.ncas-1)
    F.difProf(i,:,:,:,:)=F.pth.Profs(i+1,:,:,:,:)-F.pth.Profs(1,:,:,:,:);
end

%%
%check a few profiles

figure(1347);clf;
cas=2;
pla=8;
pin=1;
%imshow(mtomo(:,:,cas,pla),[]);

[xx,yy]=ndgrid(1:F.recsize,1:F.recsize);
    zz=zeros(F.recsize,F.recsize);
su=surf(xx,yy,zz,mtomo(:,:,cas,pla),'edgecolor','none');
colormap(gray)
angs=1:10:F.pth.n_angles;
nang=length(angs);
col=hsv(nang*4);
hold on
for pin=1:4
for ang2=1:nang;
    ang=angs(ang2);
    x=linspace(F.cen.cen(pin,1),F.Pa.ProfEndPnt(pin,ang,1),size(F.pth.Profs,5));
    y=linspace(F.cen.cen(pin,2),F.Pa.ProfEndPnt(pin,ang,2),size(F.pth.Profs,5));
    plot3(x,y,squeeze(F.pth.Profs(cas,pla,pin,ang,:)),'color',col(ang2+(pin-1)*nang,:))
    %plot3(x,y,squeeze(phanprofs(cas,pla,pin,ang,:)),'color',col(ang2+(pin-1)*nang,:))

end
end
figure(1348);clf;
dcas=1;
su=surf(xx,yy,zz,mtomo(:,:,cas,pla),'edgecolor','none');
colormap(gray)
angs=1:10:F.pth.n_angles;
nang=length(angs);
col=hsv(nang*4);
hold on
for pin=1:4
for ang2=1:nang;
    ang=angs(ang2);
    x=linspace(F.cen.cen(pin,1),F.Pa.ProfEndPnt(pin,ang,1),size(F.pth.Profs,5));
    y=linspace(F.cen.cen(pin,2),F.Pa.ProfEndPnt(pin,ang,2),size(F.pth.Profs,5));
    plot3(x,y,squeeze(F.difProf(dcas,pla,pin,ang,:)),'color',col(ang2+(pin-1)*nang,:))

end
end
%
pla=2;
figure(2356);clf;
angs=1:20:F.pth.n_angles;
nang=length(angs);
col=hsv(nang);
cas=2;
dcas=[1,2,3]
for sb=1:3
    ax=subplot(1,3,sb);
    hold on
    for ang2=2:nang-1;
        ang=angs(ang2);
        plot((1:size(F.difProf,5))*F.recres,squeeze(F.difProf(dcas(sb),pla,pin,ang,:)),...
           'color',col(ang2,:),'displayname',sprintf('%d°',ang));

    end
    if sb==1
        ylabel('\mu [1/mm]')
    end
    ylim(0.012*[0,1])
    xlabel('path length [mm]')
    title(sprintf('%3.2f kg/s',F.par.flowkgs(sb+1)))
    grid on
    plotlft=mean(sum(F.pth.Profs(sb+1,pla,pin,:,40:50)-...
        F.pth.Profs(1,pla,pin,:,40:50),5));
    disp(sprintf('mü*s film=%f',plotlft))
end
legend('location','best')

%
figure(2357);clf;
angs=1:20:F.pth.n_angles;
nang=length(angs);
col=hsv(nang);

tstr={'no flow','large flow'};
for sb=1:3
    ax=subplot(1,3,sb);
    hold on
    for ang2=2:nang-1
        ang=angs(ang2);
        plot((1:size(F.difProf,5))*F.recres,squeeze(F.pth.Profs(sb,pla,pin,ang,:)),...
           'color',col(ang2,:),'displayname',sprintf('%d°',ang));
    end
    if sb==1
        ylabel('\mu [1/mm]')
    end
    %ylim(1e-5*[-0.5,2.5]/ratio/16)
    xlabel('path length [mm]')
    title(sprintf('%3.2f kg/s',F.par.flowkgs(sb+1)))
    grid on
end
legend('location','best')



%
% 1st approach: sum up the integrals
F.LFTraw=sum(F.pth.Profs(:,:,:,:,F.pth.lftSumRange),5)*F.recres;%... % sum up the grey value
for i =1:3
    F.LFT(i,:,:,:)=(F.LFTraw(i+1,:,:,:)-F.LFTraw(1,:,:,:))/F.XS;
end

save('F','F')

% overview plot
for cas=[1,2,3]
fh=figure(5683+cas);clf;
for sb=1:4
    subplot(1,4,sb);
    im=convn(squeeze(F.LFT(cas,:,sb,:)),mt.ElliKernelGen(1,5,1,2),'same');
    imagesc(im)
    %caxis(1e-5*[0,18]*fixfac)
    caxis(max(F.LFT(:))*[0 1])
    title(sprintf('[%d %d]',cas,sb))
    
    set(gca,'YDir','normal')
    %set(gca,'XDir','reverse')
    %xticks([1 91 81]);
    yticks([]);
    xticks([1 90 181])
    xticklabels({'-90°','0°','90°'})
    
end
colorbar
end

%% make nice lft map overview plot
pub.LFTOverview(F,T)

%% make zboray-like plot
pub.LFTlikeZboray(F,T)


%% make the story publication plots
% path vizualisation
pub.paths1(F,T,mtomo)

%% make joint zboray damsohn dotox lft profile plot
pub.ZDD(T,F,C)
pub.ZDDscaled(T,F,C)
%% check or find a good F.pth.lftSumRange
pub.sumrange(F,T,mtomo)

%% make neutrons vs xray plot
pub.NeutronXray(T,F,mtomo)

%% CAD alignment
pub.CADAlignment(F,T,mtomo,F.cen.P2)

%%
load('E:\20171213 final campaign\4_rec_01_01.mat') % appears as 'd'
dfilm=f.loadSingleVariableMATFile('E:\20171213 final campaign\4_rec_06_01.mat');
%%
pub.misalignment(F,T,dfilm(:,:,400)-d(:,:,400))

%%

%% load in an empty tomo as d
pub.lftsurfplot(F,T,d,mtomo)

%% vanes plot
pub.vanesplot(T,F,d,mtomo)

%% profiles plot
pub.lftprofiles(F,T)

%% check if the symmetry is valid
pub.lftprofiles2(F,T)

%% calvin scaled plots
pub.Calvin(F,T,C)
pub.CalvinNoSpacer(F,T,C,Z)

%% 360° plot
pub.lft360profiles(F,T)

%% make pin centered plot
pub.PinCentered(F,T)

%% check if the alignment is correct
figure(1745);clf
subplot(1,2,1)
imshow(d(:,:,310)',[]);
set(gca,'YDir','normal')

subplot(1,2,2)
imshow(squeeze(mtomo(:,:,2,24)-mtomo(:,:,1,24))',[]);
set(gca,'YDir','normal')

%% some tomos for publication
corr=f.loadSingleVariableMATFile('E:\20171213 final campaign\1_corr_01_01.mat');
pub.ICONEtomos(F,T,mtomo,corr)
%%
pub.ICONErad(F,T,mtomo,corr)
%% plots for the validation section
%% plot roberts simulations for 0.1, 0.3, 0.5 mm
pub.RobertTomo(F,T,mtomo)

%% make a histogram of the data

% cas=6;rep=1;T.cas=cas;T.rep=rep;
% d=f.ReadXrayTomo(T);
pub.DataHist(d,T,F)


%% void faction
% prepare phantoms & masks
F.void.description='stuff related to the void fraction';
[F.void.phan,F.void.phanCoord]=val.makePhantomSeries(0.8,269,0.127);
% make a mask for the sub channel
F.void.subchan=(abs(F.void.phanCoord(:,:,1))+abs(F.void.phanCoord(:,:,2)))>6.7*sqrt(2);
figure(8499);clf
imshow((F.void.subchan+F.void.phan)',[]);
set(gca,'YDir','normal');
title('subchannel box')
% make a mask for the quadrants
F.void.quadrants=zeros(269);
ii=[-4:3]; % for 8-fold division
for i=1:length(ii)
    F.void.quadrants(atan2(F.void.phanCoord(:,:,2),...
        F.void.phanCoord(:,:,1))>ii(i)*pi/4)=mod(i+3,length(ii))+1;
end
% make mask plots
figure(8500);clf
imshow((F.void.quadrants+F.void.phan)',[]);set(gca,'YDir','normal');title('quadrants')
figure(8502);clf
imshow((F.void.quadrants+F.void.subchan+F.void.phan)',[]);set(gca,'YDir','normal');title('all')
F.void.finalmask=(F.void.phan==3).*(F.void.subchan==0).*F.void.quadrants;
figure(8503);clf
imshow((F.void.finalmask)',[]);set(gca,'YDir','normal');title('final')
% calculate voids
F.void.voidmeanraw=nan(length(ii),4,F.h);
F.void.mean=nan(length(ii),3,F.h);
for cas=1:4
    
    for pla=1:F.h
        tom=mtomo(:,:,cas,pla);
        for reg=1:length(ii)
            F.void.voidmeanraw(reg,cas,pla)=mean(tom(F.void.finalmask==reg));
        end
    end
    if cas>1
        F.void.mean(:,cas-1,:)=F.void.voidmeanraw(:,cas,:)-F.void.voidmeanraw(:,1,:);
    end
end
F.void.mean(F.void.mean<=0)=nan
clear tom


%% plot void profiles
pub.voidprofile(F,T,mtomo)


%% control plot
figure(2378);clf;
cas=2;
pla=34;
imshow(((mtomo(:,:,cas,pla)-mtomo(:,:,1,pla)).*(F.void.finalmask~=0))',[]); set(gca,'YDir','normal');

%% marke radiograph for thesis
T.cas=1;
T.rep=1;
d=f.ReadXrayTomo(T);
%%
pub.radiograph(T,F,d(:,:,50))


%% make time signal for one row
% read in a radiograph video and do stuff then
if 1
T.cas=2;
T.rep=1;
T.d.List{2,1}='20171213_rad_06_45_1.seq';
T.d.nframes(T.cas,T.rep)=300;
drad=f.ReadXrayTomo(T);
T.d.List{2,1}=[];
end
%% integral void
F.ivoid.mask=val.make_phantom(0.5,269,0.127);
F.ivoid.mask(F.ivoid.mask<3)=0;
F.ivoid.mask(F.ivoid.mask>=3)=1;
%control plot
figure(26780);clf;imshow(F.ivoid.mask',[]);set(gca,'ydir','normal')
figure(26781);clf;imshow(cat(3,F.ivoid.mask,...
    f.imnorm(mtomo(:,:,2,30)-mtomo(:,:,1,30)),zeros(269)),[]);set(gca,'ydir','normal')

F.ivoid.mask=pub.ImSegExp(F.ivoid.mask,4);
figure(26782);clf;imshow(cat(3,F.ivoid.mask,...
    f.imnorm(mtomo(:,:,2,30)-mtomo(:,:,1,30)),zeros(269)),[]);set(gca,'ydir','normal')

F.ivoid.description='integral void / void fractino of the full channel cross section'
F.ivoid.ivoid=nan(3,63);
for cas=1:3
    for pla=1:F.h
        tom=mtomo(:,:,cas+1,pla)-mtomo(:,:,1,pla);
        F.ivoid.ivoid(cas,pla)=mean(tom(F.ivoid.mask==1))/F.XS;
        

    end
end
%% make the figure
pub.ZborayVoid(T,F)


%% tomography example for presentation
pub.TomExample(T,F,mtomo)

%% try the center of mass of the void %% didnt work or doesnt show interesting stuff
F.void.voidcenraw=nan(length(ii),4,F.h,2);
F.void.voidcen=nan(length(ii),3,F.h,2);
[xx,yy] = ndgrid(F.p.xcoord,F.p.xcoord);
for cas=1:4
    for pla=1:F.h
        tom=mtomo(:,:,cas,pla);
        for reg=[1 3 5 7]

            mask=(F.void.finalmask==reg)+(F.void.finalmask==reg+1);
            F.void.voidcenraw(reg,cas,pla,1)=...
                sum(mask(:).*tom(:).*xx(:))/sum((mask(:).*tom(:)));
            F.void.voidcenraw(reg,cas,pla,2)=...
                sum(mask(:).*tom(:).*yy(:))/sum((mask(:).*tom(:)));
        end
    end
    if cas>1
        F.void.voidcen(:,cas-1,:,:)=F.void.voidcenraw(:,cas,:,:)-F.void.voidmeanraw(:,1,:,:);
    end
end
%%
figure(256368);clf
imshow(tom',F.p.RA',[]);
set(gca,'YDir','normal')
hold on

col=jet(F.h);
cas=1;
for reg=[1 3 5 7]
for pla=6:F.h-1
    p(pla)=plot(squeeze(F.void.voidcen(reg,cas,pla,1)),...
         squeeze(F.void.voidcen(reg,cas,pla,2)),...
        'x','color',col(pla,:),'displayname',sprintf('pla %d',pla));
end
end
grid on
legend(p(22:10:F.h))

%% make a video with spacer
corr=exp(-f.loadSingleVariableMATFile('E:\20171213 final campaign\1_corr_06_01.mat'));
%%
pub.spcerVid(F,T,corr)


%%
figure(1346);clf
imshow((mtomo(:,:,cas,pla)-mtomo(:,:,1,pla))'+1e-5*vmask,[])

void=(mtomo(:,:,cas,pla)-mtomo(:,:,1,pla))'.*vmask;
figure(257);clf;
imshow(void,[])
title('void signal')

%% 
figure(4835);clf
ax=cla();
hold on
 cases=[2 3 4]
tstr={'high flow','mid flow','low flow'} 
for i=1:3 
subplot(1,3,i)
hold on
for reg=inner
    plot(squeeze(mean(reg,cases(i),nospac)-mean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end

grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title(tstr{i})
ylabel('entrainment')
end
%% de-centering plot


cas=8;
rep=1;
fname=sprintf('%s%02d_%02d.mat','1_corr_',cas,rep);
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


%%
pla=600;
[ax,fh]=pub.decenter(T,d);

%% center shift  scan

sino=squeeze(d(:,pla,:));

% center exact
%baseshift(plane)
nshift=41; % number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);


%%
mask=f.MakeTomoMask('fan',f.fraccircshift(...
                sino,extrashift(21)),T.Rec.recsize,deg2rad(T.Rec.angles),...
                T.Rec.detPitch,T.Rec.src,T.Rec.det,...
                T.Cen.MaskThresh,0.000005,[150 150]);
            figure(23467);clf;imshow(mask,[])
            
 %%           
% reconstruct
rec=zeros(T.Rec.recsize,T.Rec.recsize,nshift);
for i =1:nshift
rec(:,:,i)=a.FBPexplFan(...
    f.fraccircshift(sino,extrashift(i))',...
    T.Rec.recsize,deg2rad(T.Rec.angles),T.Rec.detPitch,T.Rec.src,T.Rec.det);
    [q(i),~]=qc.VarianceQuality(rec(:,:,i),...
        mask);
end
%%
figure(23451);clf
imshow(mask,[])
title('mask')
fh=figure(2314);clf
mask2=mask;
mask2(round(T.Rec.recsize/2):end,1:round(T.Rec.recsize/2))=0;
imshowpair(imgradient(rec(:,:,21))',mask2'>1);
set(gca,'YDir','normal')
text(153,250,'image gradient','Color',[0 1 0],'FontSize',12)
text(220,230,'mask','Color',[1 0.2 1],'FontSize',12)
ax=gca;
ax.Units='pixel';
ax.Position=[0 0 T.Rec.recsize T.Rec.recsize];
fh.Units='pixel';
fh.Position=[50 50 T.Rec.recsize T.Rec.recsize];
fh.Units='centimeters';
fh.PaperPosition=fh.Position;
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperUnits', 'centimeters')
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

fname=sprintf('Qualitymask1');
print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
%open(sprintf('%s%s.pdf',T.F.saveto,fname))
%%
fh=figure(4894);clf;
plot(extrashift,q,'x')
grid on
ylh=xlabel('shift [pixel]');
xlh=ylabel('$\sigma(\nabla Im\mid_{Mask})$','Interpreter','latex');
pub.BFxlab(xlh,T.F)
pub.BFylab(ylh,T.F)


fh.Units='centimeters';
fh.Position=[2 2 7.1173 7.1173];

set(fh, 'PaperUnits', 'centimeters')
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

ax = gca;
pub.BFaxis(ax,T.F)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fname=sprintf('Qualitymask2');
print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.F.saveto,fname))
%%
fh=figure(48941);clf;
plot( T.Cen.eshift(:,cas,rep))
hold on
plot( pla,T.Cen.eshift(pla,cas,rep),'x')

grid on
ylh=xlabel('z-plane');
xlh=ylabel('best shift [pixel]');
xlim([0 1024])
ylim(1.2*[-1,1])
pub.BFxlab(xlh,T.F)
pub.BFylab(ylh,T.F)


fh.Units='centimeters';
fh.Position=[2 2 7.1173 7.1173];

set(fh, 'PaperUnits', 'centimeters')
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

ax = gca;
pub.BFaxis(ax,T.F)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];
fname=sprintf('Qualitymask3');
print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.F.saveto,fname))

%%
fh=figure(49846);clf;
im=[];
xrange=30:80;
yrange=50:120;
c1=[];
for i=1:5:nshift
   im=cat(1,im,rec(xrange,yrange,i)); 
   c1=cat(1,c1,i);
end

im2=[];
c2=[];
for i=17:25
   im2=cat(1,im2,rec(xrange,yrange,i)); 
   c2=cat(1,c2,i);
end

im=cat(2,im2,im);
imshow(im',[])
set(gca,'YDir','normal')

for i=1:9
    text(25+(i-1)*51,95,sprintf('%2.1f',extrashift(c1(i))),...
        'color',[0.99 1 1],'FontSize',12,'FontName', 'Helvetica')
    text(25+(i-1)*51,25,sprintf('%2.1f',extrashift(c2(i))),...
        'color',[.99 1 1],'FontSize',12,'FontName', 'Helvetica')
end

ax=gca;
ax.Units='pixel';
ax.Position(1:2)=[0 0];

fh.Units='pixel';
fh.Position=[50 50 ax.Position(3:4)];
fh.Units='centimeters';


set(fh, 'PaperUnits', 'centimeters')
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
fh.PaperSize=fh.Position(3:4);

    fname=sprintf('Qualitymask4');
print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.F.saveto,fname))


% fh.Units='centimeters';
% fh.Position=[2 2 7.1173 7.1173];
% 
% set(fh, 'PaperUnits', 'centimeters')
% fh.PaperSize=fh.Position(3:4);
% set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
% 
% ax = gca;
% pub.BFaxis(ax,T.F)
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% fname=sprintf('Qualitymask4');
% print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
% open(sprintf('%s%s.pdf',T.F.saveto,fname))

%% old stuff




%%
load('F')
load('T')
load('recon')
d=f.loadSingleVariableMATFile('E:\20171213 final campaign\2_add_02.mat');
%% lft
lft.ProfileAnalyser(F,recon);
%%
lft.LFTMap(F)
%%
lft.plotProfiles(F)
%%
lft.plotMinMaxProfiles(F)
%%
lft.Rayplot(F,recon)
%%
lft.circleplot(F,recon)
%%
lft.plotMinMaxProfiles(F)

%%
savefast('F','F');

%% publication graphs
pub.doubletomo(F,recon)

%%
pub.axprofile(F)
%%
pub.quadtomo(F,recon)
%%
pub.horzprof(F)
%%
pub.lftmap(F)

%%
pub.quadtomo2x2(F,recon)

%%
pub.dh(F,d)
%%
pub.lftmap(F)
%%
perm=[3,1,2];
dperm=permute(d(:,31:1015,:),perm);
%%
pub.volumePlot(dperm,T)
%%
pub.printVector





%%
% new recon images
















