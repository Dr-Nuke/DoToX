% test the cone beam algo

% get T
if exist('T','var')
    clearvars -except T
    
else
    clear
    load('T')
end
clc
close all
format compact

% add astra paths
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')

% case & rep
c=1;
r=1;
d=f.BoLoad(T.fnames.corr{c,r},T);

%%
ang=deg2rad(T.Rec.angles); % projection angles list

xrange=T.Cen.start:T.Cen.stop; % rough cropping around channel
p=512;
hs=5;
yrange=p-hs:p+hs;
sinorange=T.q360.startframe:T.q360.startframe+T.q360.NFGuess-1; % 360°crop

draw=-log(d(xrange,yrange,sinorange));
shift=T.Cen.fitshift(c,p);
quadrot=70; % rotate to have the channel quadrants aligned

dwork=permute(f.fraccircshift(circshift(draw,quadrot,3),-(shift)),[1 3 2]);

detPitch=0.127; % detector pixel spacing
src=950/detPitch; % source-axis distance
det=50/detPitch; % detector-axis distance

figure(1);clf;
imshow(a.FBPexplFan(dwork(:,:,5)',319,ang,0.127,950,50),[])
figure(2);clf;
imshow(a.FBPexpl(dwork(:,:,5)',319,ang),[])


%%

volxpix=50;
volypix=50;
volzpi=70;
det_pitch=0.127; % detector pixel pitch
npix=319 % number of detector pixels

vol_geom=astra_create_vol_geom(volxpix,volypix,volzpix);
rec_id=astra_mex_data3d('create', '-vol', vol_geom);
proj_geom = astra_create_proj_geom('fanfalt', 1.0, 1.0, 128, 192, angles);




rec=astre_mex_data3d('get',rec_id);
astra_mex_data3d('delete',rec_id);

astra_mex_data3d('info');


