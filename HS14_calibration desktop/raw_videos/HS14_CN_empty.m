clear all 
clc
format compact
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

tic;
fdir='C:\data\Tomo_HS14\02_rawdata\empty\'

nmax=406;


nx=1200;
ny=2450;

xpix=1:4:nx;
ypix=1:4:ny;

data=zeros(length(xpix),length(ypix),nmax);
FNamStr={'empty_','.fits'};


for i=1:nmax
    disp(i)
    FNum=sprintf('%04u',i);
    
    fname=strcat(FNamStr{1},FNum,FNamStr{2});
    fpath=strcat(fdir,fname);
    %disp(fpath)
    a=fitsread(fpath)';
    a=a(xpix,ypix);
    data(:,:,i)=a;
    
    %fig=figure(4);
    %imshow(a,);set(gca,'YDir','normal')

    
    
end

contrast=[40000 58000]
%%

xdose=1:80;
ydose=1:length(ypix);
meandose=mean(mean(mean(data(xdose,ydose,:))));
checkdose=mean(mean(data(xdose,ydose,:)));
data2=data(:,:,checkdose>48000);


%% video


for i=1:size(data2,3)
    disp(i)
    fig=figure(4);
    imshow(squeeze(data2(:,:,i))',contrast);set(gca,'YDir','normal')    

    drawnow()
    M(i)=getframe(fig);
    
end

vidpath='G:\cbolesch\Data\';
vidname='HS14 CN empty';

video=VideoWriter([strcat(vidpath,vidname)],'MPEG-4');
video.FrameRate=25;
open(video)
writeVideo(video,M);
close(video)


t2=toc
disp(sprintf('run time: %.0f seconds',t2))