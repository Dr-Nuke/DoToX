clear all 
clc
format compact
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

tic;
fdir1='G:\cbolesch\Data\P20150733 Dezember 2015 tomo\02_rawdata\06_transiente\';


nmax1=45;
nmax2=1514;

nx=1150;
ny=2560;

data=zeros(288,640,nmax1+nmax2);
FNamStr1={'transiente__','.fits'};
FNamStr2={'p40_T90_2_','.fits'};


for i=1:nmax1
    f_BoCount(i,20,15,5);
    FNum=sprintf('%04u',i);
    
    fname1=strcat(FNamStr1{1},FNum,FNamStr1{2});
    fpath1=strcat(fdir1,fname1);
    %disp(fpath)
    a=fitsread(fpath1)';
    a=a(1:4:nx,1:4:ny);
    data(:,:,i)=a;
    
    %fig=figure(4);
    %imshow(a,[0,30000]);set(gca,'YDir','normal')
end
% %%
% for i=1:nmax2
%     f_BoCount(i+nmax1,20,15,5);
%     FNum=sprintf('%04u',i);
%     
%     fname2=strcat(FNamStr2{1},FNum,FNamStr2{2});
%     fpath2=strcat(fdir2,fname2);
%     %disp(fpath)
%     a=fitsread(fpath2)';
%     a=a(1:4:nx,1:4:ny);
%     data(:,:,i+nmax1)=a;
%     
%     %fig=figure(4);
%     %imshow(a,[0,30000]);set(gca,'YDir','normal')
% end


%%

xdose=20:70;
ydose=30:600;
meandose=mean(mean(mean(data(xdose,ydose,:))));
checkdose=mean(mean(data(xdose,ydose,:)));
data2=data(:,:,checkdose>20000);


%% video


for i=1:size(data2,3)
    f_BoCount(i,20,15,5);
    
    fig=figure(4);
    imshow(squeeze(data2(:,:,i))',[0,30000]);set(gca,'YDir','normal')
    M(i)=getframe(fig);
    
end
viddir='G:\cbolesch\Data\';
vidname='H515 ICON transiente';

video=VideoWriter([strcat(viddir,vidname)],'MPEG-4');
video.FrameRate=5;
open(video)
writeVideo(video,M);
close(video)


t2=toc
disp(sprintf('run time: %.0f seconds',t2))