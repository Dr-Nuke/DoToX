% minimal working example for phantom-forwardprojection-reconstruction
% via astra fanbeam & parallel beam
%clear

if 0  % load T & d
    clear
    close all
    cas=1;rep=1;T.cas=cas;T.rep=rep;
    load('T')
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
    
    v=nan(niter,BS(2));
    xlist=nan(niter,BS(2));
    eshift=nan(1,BS(2));
    recon=zeros(recsize,recsize,BS(2),'single');
end
%%
if 1
    
    [phan,coord]=val.make_phantom(0,3000,0.0127);
    phan=circshift(phan,-10,1); % x shift
    phan=circshift(phan,-15,2); % y shift
    phan(phan==0)=0.7; % al
    phan(phan==1)=0; %air
    phan(phan==2)=0.18; %h20
    phan(phan==3)=0.7*4.475/1411; %cl3 vap
    phan_emp=phan;
    phan(phan==4)=0.78; %cl3
    phan=phan/10;
    phan=imresize(phan,[300,300]);
    phanrange=16:284;
    phan=phan(phanrange,phanrange);
    
    phan_emp(phan_emp==4)=0.78*4.475/1411; %cl3 vap
    phan_emp=imresize(phan_emp,[300,300]);
    phan_emp=phan_emp/10;
    phan_emp=phan_emp(phanrange,phanrange);
end
%
% parameters
% d_ph=30; % phantom area width in mm
% n_ph=269; % Phantom area width in pixel
% v_ph=1;   % phantom pixel value in 1/cm
% v_ph=v_ph/10; % 1/cm to 1/mm
% p_ph=d_ph/n_ph; % phantom pixel pitch in mm

d_ph=30; % phantom area width in mm
n_ph=269; % Phantom area width in pixel
v_ph=1;   % phantom pixel value in 1/cm
v_ph=v_ph/10; % 1/cm to 1/mm
p_ph=0.127; % phantom pixel pitch in mm
d_ph=p_ph*n_ph; % phantom area width in mm

n_d=269;  % number of detector pixels
p_d=0.127;   % detector pixel spacing in mm

p_r=0.127;    % reconstruction pixel size in mm
d_r=30;            % reconstruction area width in mm
n_r=round(d_r/p_r); % reconstruction area width in pixel
n_r=269;

n_ang=1357; % number of projection angles / Question 1 / see line 113
src=950;  % source-origin-distance in mm
dets= 50;   % detector-origin distance in mm
det=50;  %dets(i);
src=1000-det;

%phantom
V = zeros(n_ph);
xrang=round(n_ph/4):round(n_ph/4*3);
yrang=xrang;
V(xrang,yrang)=v_ph; % set some phantom part
V=phan_emp; % 
%RV=phanweighted(:,:,1);

%angles
ang=linspace(0,2*pi,n_ang+1); % make the angle list in 2 steps
ang(end)=[];

mode={'parallel','fanflat','real (fan)'}; % for plot titles

% phantom plot
figure(1) ;clf
imagesc(V');set(gca,'YDir','normal')
title(sprintf('phantom, value=%2.1f/mm or %2.1f/cm',v_ph,10*v_ph))
colorbar

% creat sinogram. use provided example code, slightly modiefied
clear sinogram sino2
sinogram(:,:,1)=val.parforward_new(V,n_ph,d_ph,p_d,n_d,ang);
sinogram(:,:,2)=val.fanforward_new(V,d_ph,p_d,n_d,ang,src,det);
sinogram(:,:,3)=squeeze(d(:,400,:))';

%sinogram plot
figure(2) ; clf
for i =1:size(sinogram,3)
    ax(i)=subplot(1,size(sinogram,3),i);
    imshow(sinogram(:,:,i), []);set(gca,'YDir','normal')
    colorbar
    tstr=strcat('$\int \mu \cdot s ds $ ',sprintf(' %s',mode{i}));
    title(tstr,'Interpreter','latex')
    xlabel(sprintf('straight path integral %f',max(sinogram(1,:,i))))
    sino2(:,:,i)=exp(-sinogram(:,:,i));
    xlabel(sprintf('%f',max(sinogram(1,:,i))))
end
linkaxes(ax);

figure(4) ; clf
for i =1:size(sinogram,3)
    subplot(1,size(sinogram,3),i)
    sino2(:,:,i)=exp(-sinogram(:,:,i));
    imshow(sino2(:,:,i), []);set(gca,'YDir','normal')
    colorbar
    tstr=strcat('$e^{-\int \mu \cdot s ds} $ ',sprintf(' %s',mode{i}));
    title(tstr,'Interpreter','latex')
    %xlabel(sprintf('straight path integral %f',max(sinogram(1,:,i))))
    
end
%
% creat reconstruction. use provided example code, slightly modiefied
clear rec
rec(:,:,1)=val.par_new(sinogram(:,:,1),ang,n_r,p_r,p_d);
rec(:,:,2)=val.fan_new(sinogram(:,:,2),ang,src,det,n_r,p_r,p_d);
rec(:,:,3)=val.fan_new(sinogram(:,:,3),ang,src,det,n_r,p_r,p_d);

% recon plot
figure(3) ;clf;
for i =1:size(rec,3)
    subplot(1,size(rec,3),i)
    % find the mean of a square that represents the phantom's unity region
    ev1=max(min(round((p_ph/p_r*(-n_ph/2+xrang(5))+n_r/2)),n_r),0);
    ev2=max(min(round((p_ph/p_r*(-n_ph/2+xrang(end-5))+n_r/2)),n_r),0);
    evarang=ev1:ev2;
    tf(i)=mean(mean((rec(evarang,evarang,i))));% "transfer function"
    imagesc(rec(:,:,i)',[0,0.1]);set(gca,'YDir','normal')
    hold on
    colorbar
    rectangle('Position',[ev1,ev1,ev2-ev1+1,ev2-ev1+1])
    title(sprintf('scaled recon %s;\n mean of box =%f',mode{i},tf(i)))
end

% sinogram proifiles
figure(5);clf
gca()
hold on
for i=1:size(sinogram,3)
    plot(sinogram(1,:,i))
end
legend(mode)
title('sinogram profiles')
grid on

%econ profile plot
figure(52);clf
gca()
hold on
for i=1:size(rec,3)
    plot(rec(139,:,i))
end
legend(mode)
title('recon profiles')
grid on

figure(6) %imdif plot
imshowpair(rec(:,:,2)',rec(:,:,3)');set(gca,'YDir','normal')
title('fan phantom vs. real')
figure(7) %imdif plot
imshowpair(rec(:,:,2)',mtomo(:,:,1,18));set(gca,'YDir','normal')
title('fan phantom vs. Micha')

%%
% dtstance evaluation plot
ev1=max(min(round((p_ph/p_r*(-n_ph/2+xrang(5))+n_r/2)),n_r),0);
ev2=max(min(round((p_ph/p_r*(-n_ph/2+xrang(end-5))+n_r/2)),n_r),0);
evarang=ev1:ev2;
dets=[0 10 20 50 100 200 300 400 500];
tf=zeros(1,length(dets));
for i=1:length(dets)
    det=dets(i);
    src=1000-det;
    sinogram(:,:,i)=val.fanforward_new(V,d_ph,p_d,n_d,ang,src,det);
    rec(:,:,i)=val.fan_new(sinogram(:,:,i),ang,src,det,n_r,p_r,p_d)...
        ;
    tf(i)=mean(mean((rec(evarang,evarang,i))));% "transfer function"
end
%
figure(165979);clf
plot(dets,tf/v_ph,'x-')
grid on
hold on
plot(dets(4),tf(4)/v_ph,'ro');
%
annotation('textarrow',[0.4,0.25],0.85*[0.9 0.98],'String',sprintf('real geometry: %f',tf(4)))
xlabel('detector-origin')
title(sprintf('variation at constant src-det=%dmm',1000))
ylabel('transfer function')
