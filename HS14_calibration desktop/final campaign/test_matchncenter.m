% here we try to match 
% get data
clear 
close all
load('T')
d1=f.loadSingleVariableMATFile('E:\20171213 final campaign\1_corr_01_01.mat');
d8=f.loadSingleVariableMATFile('E:\20171213 final campaign\1_corr_08_01.mat');

d1=d1(173:491,:,:);
d8=d8(173:491,:,:);

%%
close all
plane = 100;
height=1;
thresh=0.002;
ang=T.Rec.angles;
% log, squeez, and mean of stack height
im1=-log(squeeze(mean(d1(:,plane:plane+height-1,80:1357+80-1),2)));
im8=-log(squeeze(mean(d8(:,plane:plane+height-1,80:1357+80-1),2)));

% shift to have the quadrants aligned
im12=circshift(im1,T.Rec.SinoRotShift,2);
im82=circshift(im8,T.Rec.SinoRotShift,2);



%create the unique mask
im13=f.fraccircshift(im12,-(T.Cen.fitshift(1,100)));
rec1=iradon(im13,ang,'spline','Hann',1,319);
[~,mask,gm1]=qc.kern(rec1,thresh,[160,160]);

figure(15);clf;
imshow(cat(1,rec1/max(rec1(:)),gm1/max(gm1(:)),mask),[])
title('recon - grad - mask')



%% check centering
nshift=12; % number of tests
minshift=-0.5; % in pixels
maxshift=1.5; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);

plane=640;
xrange=97:197;
yrange=20:120;

% log, squeez, and mean of stack height
im1=-log(squeeze(mean(d1(:,plane:plane+height-1,80:1357+80-1),2)));
im8=-log(squeeze(mean(d8(:,plane:plane+height-1,80:1357+80-1),2)));

% shift to have the quadrants aligned
im12=circshift(im1,T.Rec.SinoRotShift,2);
im82=circshift(im8,T.Rec.SinoRotShift,2);



totim=[];
for i=1:nshift
% centering shift
im13=f.fraccircshift(im12,-(T.Cen.fitshift(1,100)+extrashift(i)));
im83=f.fraccircshift(im82,-(T.Cen.fitshift(4,100)+extrashift(i)));

%reconstruction
rec1(:,:,i)=iradon(im13,ang,'spline','Hann',1,319);
rec8(:,:,i)=iradon(im83,ang,'spline','Hann',1,319);

[v(i,1),g1(:,:,i)]=qc.VarianceQuality(rec1(:,:,i),mask);
[v(i,2),g8(:,:,i)]=qc.VarianceQuality(rec8(:,:,i),mask);


end
b=[];
d=[];
for i=1:4
    a=[];
    c=[];
    for j=1:3
        k=j+(i-1)*3;
        a=cat(2,a,g8(:,:,k));
        c=cat(2,c,rec8(:,:,k));
    end
    a=cat(2,a,zeros(319,319));
    b=cat(1,b,a);
    c=cat(2,c,zeros(319,319));
    d=cat(1,d,c);
end


figure(1);clf;
imshow(b,[])
title('centering opt.')

figure(10);clf;
imshow(d,[])
title('centering opt.')

figure(2);clf;
plot(extrashift,v(:,1),'xr','displayname','no film' )
hold on
plot(extrashift,v(:,2),'xb','displayname','withfilm')
xlabel('centering shift')
ylabel('gradient mask size')
grid on
title('optimization for centering')
legend()

% csvwrite('WOFilm10SliceAverage.csv',rec1);
% csvwrite('WFilm10SliceAverage.csv',rec8);

%%




% check angles

nang=12;
minang=-1;
maxang=+1;
extraang=linspace(minang,maxang,nang);

im15=f.fraccircshift(im12,-(T.Cen.fitshift(1,100)));
im85=f.fraccircshift(im82,-(T.Cen.fitshift(4,100)));
totimang=[];
for i = 1:nang
    ang=linspace(0,360+extraang(i),T.q360.nFrames(1,1)+1);
    ang(end)=[];
    
    rec15(:,:,i)=iradon(im13,ang,'spline','Hann',1,319);
    rec85(:,:,i)=iradon(im83,ang,'spline','Hann',1,319);

[vang(i,1),gang1(:,:,i)]=qc.VarianceQuality(rec15(:,:,i),mask);
[vang(i,2),gang8(:,:,i)]=qc.VarianceQuality(rec85(:,:,i),mask);


end

%make image
b=[];
for i=1:4
    a=[];
    for j=1:3
        k=j+(i-1)*3;
        a=cat(2,a,gang1(:,:,k));
    end
    a=cat(2,a,zeros(319,319));
    b=cat(1,b,a);
end


figure(3);clf;
imshow(b,[])
title('angle opt.')

figure(4);clf;
plot(extraang,vang(:,1),'xr','displayname','no film' )
hold on
plot(extraang,vang(:,2),'xb','displayname','withfilm')
xlabel('extra angle')
ylabel('gradient mask size')
title('optimization for angle')
grid on
legend()


%%
nshifft=21;
smin=-2;
smax=2;





























