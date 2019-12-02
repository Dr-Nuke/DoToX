clear all 
clc
format compact
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

tic;
fdir='G:\cbolesch\20151208 CNHT X-ray campaign\'
fname='5 boiling after CHCL3 pump trip.avi'

%%
v = VideoReader(strcat(fdir,fname));
mov=read(v);
mov2=mov(:,:,1,:);

vidname='bubbly'

video=VideoWriter([strcat(fdir,vidname)],'MPEG-4');
video.FrameRate=25;
open(video)
writeVideo(video,mov2);
close(video)