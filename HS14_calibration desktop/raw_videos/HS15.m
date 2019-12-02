clear all 
clc
format compact
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

tic;
fdir='C:\data\P20150733 Dezember 2015\02_rawdata\05_p40_T90_2\'

nmax=1514;


nx=1150;
ny=2560;

data=zeros(288,640,nmax);
FNamStr={'P40_T90_2_','.fits'};


for i=1:nmax
    
    FNum=sprintf('%04u',i);
    
    fname=strcat(FNamStr{1},FNum,FNamStr{2});
    fpath=strcat(fdir,fname);
    %disp(fpath)
    a=fitsread(fpath)';
    a=a(1:4:nx,1:4:ny);
    data(:,:,i)=a;
    
    fig=figure(4);
    imshow(a,[0,30000]);set(gca,'YDir','normal')

    
    
end
    
%%

xdose=20:70;
ydose=30:600;
meandose=mean(mean(mean(data(xdose,ydose,:))));
checkdose=mean(mean(data(xdose,ydose,:)));
data2=data(:,:,checkdose>22500);


%% video


for i=1:size(data2,3)
    disp(i)
    
    fig=figure(4);
    imshow(squeeze(data2(:,:,i))',[0,30000]);set(gca,'YDir','normal')
    M(i)=getframe(fig);
    
end

vidname='chcl3'

video=VideoWriter([strcat(fdir,vidname)],'MPEG-4');
video.FrameRate=25;
open(video)
writeVideo(video,M);
close(video)


t2=toc
disp(sprintf('run time: %.0f seconds',t2))