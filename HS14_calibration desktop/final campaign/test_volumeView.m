d1=f.loadSingleVariableMATFile('E:\20171213 final campaign\4_rec_01_01.mat');
d2=f.loadSingleVariableMATFile('E:\20171213 final campaign\4_rec_06_01.mat');
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/apriori/')
P=f42_ChannelPhanTransform()
%%
im1=d1(:,:,350:end);
im2=d2(:,:,350:end);

dim=im2-im1;
fh=figure(1231);clf 


imshow(imgradient(squeeze(d2(:,:,100)))',[])
set(gca,'Ydir','normal')
hold on

colorbar
axis on



P2=f42_RotPhan(P,deg2rad(45.5));
P2=f42_PhanZoom(P2,7.9);
P2=f42_PhanMove(P2,134.8,134.4);
f42_PlotPhan(P2,gca,[1 0 0],[1 4])

set(fh,'units','normalized',...
    'outerposition',[0 0 1 1])

