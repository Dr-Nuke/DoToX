%% Tabula Rasa
clear all; 
clc; 
format compact; 
close all;

%% Open Files
%fid=fopen('E:\ICONchris_Xray\Xray\Tuning\July2nd\Flatpanel\Dynamics85fps2x2bin\slowspeed.seq');
filepath=('C:\Users\robersl\Desktop\Working Folder\110kV 8mA 10fps 3mmAl 0_4mmCu ensi tomo.seq');
fileInfo = dir(filepath);
imsize=[1024 640];
header_size=2048;
fileSize = fileInfo.bytes;
fileSize=floor((fileSize-header_size)/imsize(1)/imsize(2)/2);
fid=fopen(filepath);
hed = fread(fid,header_size);%the header size migth need to be adjusted depending on image settings
%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
img=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16')));
fclose(fid);

% Create a 3D Matrix out of the imported data
img=reshape(img,[1024 640 fileSize]);
img_original=img;
img=double(img);

%% Dose Correction
img_dose_correction=img_original;
img_dose_correction(1:50,:)=NaN;
img_dose_correction(:,230:485)=NaN;
img_dose_correction(345:810,1:210)=NaN;
img_dose_correction=double(img_dose_correction);
img_dose_correction=squeeze(nanmean(nanmean(img_dose_correction,1),2));

for i=1:size(img,3)
    img(:,:,i) = img(:,:,i)./img_dose_correction(i);
end

%% Calculate the Flat field
% Cut Picture to region of interest (i.e. Remove frames at beginning and the end and cut channel
tic
img_ref=zeros(size(img,1),size(img,2));

img_fit=img;
img_fit=img_fit(:,:,20:end-20);
img_fit=mean(img_fit,3);
img_fit(1:50,:)=NaN;
img_fit(:,230:485)=NaN;
img_fit(345:810,1:210)=NaN;

% Plot Mean Background to see if it's good
figure(1)
imagesc(img_fit)
caxis([0 1.5])
colorbar
%%
% Create point (x,y,z) of the background data for the 2D-fit
x=zeros(numel(img_fit),1);
y=zeros(numel(img_fit),1);
z=zeros(numel(img_fit),1);
k=1;
% Disregard edges because of errors
for i=5:size(img_fit,1)-5
    for j=2:size(img_fit,2)-1
        if img_fit(i,j) > 1
            x(k)=i;
            y(k)=j;
            z(k)=img_fit(i,j);
            k=k+1;
        end
    end
end
x=x(1:k-1);
y=y(1:k-1);
z=z(1:k-1);
figure(2)						% Check if data looks good
plot(z)
f = fit( [x,y],z, 'poly22' );	% Do the fit
figure(3)						% Show the fit
plot(f,[x(1:100:end),y(1:100:end)],z(1:100:end))

% Create Reference image from fit
for i=1:size(img_ref,1)
    for j=1:size(img_ref,2)
        img_ref(i,j)=f(i,j);
    end
end
%%
% Normalize each frame with the Background
for i=1:size(img,3)
    img(:,:,i) = img(:,:,i)./img_ref;
end

%% Reduce image to the region of interest
cropx=270:1024-10;
cropy=220:490;
img=img(cropx,cropy,:);

%fid=fopen('G:\cbolesch\20151208 CNHT X-ray campaign\78 boiling p00 h1_0000 h2_ 0000 30fps rot.seq');
%hed = fread(fid,2048);
%imgref=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16')));
%imgref=reshape(imgref,[1024 640 fileSize]);

%% 250-480
plane=140; %row to reconstruct
for i = 1:230	% For a Video go through different rows
   plane=i+249;
steps=903;		% Number of projections
startstep=47;	% Choose starting frame
dtheta=360/(steps-1);	% Angle between projections
%dtheta = 3.98;
%j=99;%to start the fil sino for 76 boiling p40
%j=120;%to start the fil sino for 75 boiling p40
%for j=95:103

%for i=0:20
%dtheta = 3.91+i*0.02;

% Create Sinogram from data 
sino=double(squeeze(img(plane,:,:)));
sino=sino(:,startstep:(startstep+steps-1));

%for 77 boiling p40
%sino=sino(:,1:i);sinoref=sinoref(:,j:(j+i-1));dtheta=180/(i-1);%for 76 boiling p40
%sino=sino(:,j:(j+i-1));sinoref=sinoref(:,1:i);dtheta=180/(i-1);%for 75 boiling p40

%Dose normalisation
% Di=mean(sino(end-20:end,:),1);
% Dref=mean(sinoref(end-20:end,:),1);
% Di=repmat(Di,size(sino,1),1);
% Dref=repmat(Dref,size(sino,1),1);
% pro2=-log(Dref.*sino./Di./sinoref);
%pro2=-log(sinoref./Dref);
pro2=-log(sino);

% Center the projection
[c,lags]=xcov(pro2(:,1),flipud(pro2(:,end)),500,'coeff');
[~,ind]=max(c);%
beta=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/(log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));

if fix(beta)>0
    inte=fix(beta); 
    frac=beta-fix(beta);
    pro2=[pro2; repmat(pro2(end,:),abs(inte),1)]; %integer shift to center
    if mod(size(pro2,1),2) 
       [XI,YI]=meshgrid(1:size(pro2,2),1:size(pro2,1));
       pro2=interp2(XI,YI,pro2,XI,YI+frac,'cubic',mean(pro2(end,:))); %subgrid shift
    else
        [XI,YI]=meshgrid(1:size(pro2,2),1:size(pro2,1)+1);
        pro2=interp2(XI,YI,[pro2; pro2(end,:)],XI,YI+frac-0.5,'cubic',mean(pro2(end,:)));
    end %this is to get always odd sized projections, that preserves centeredness in iradon
else
    inte=fix(beta); 
    frac=beta-fix(beta);
    pro2=[repmat(pro2(1,:),abs(inte),1);pro2]; %integer shift to center
    if mod(size(pro2,1),2) [XI,YI]=meshgrid(1:size(pro2,2),1:size(pro2,1));
       pro2=interp2(XI,YI,pro2,XI,YI+frac,'cubic',mean(pro2(1,:))); %subgrid shift
    else [XI,YI]=meshgrid(1:size(pro2,2),1:size(pro2,1)+1);
     pro2=interp2(XI,YI,[pro2; pro2(end,:)],XI,YI+frac-0.5,'cubic',mean(pro2(1,:)));
    end %this is to get always odd sized projections, that preserves centeredness in iradon
end

%figure();imagesc(pro2)
%
fig=figure(4);
axis equal
set(fig, 'Position', [2148 188 738 420])
% Calculate the tomography
I=iradon(pro2(:,1:steps),dtheta,'spline','Han',250); 
%imbo4(I)
subplot(1,3,[1:2])
imagesc(squeeze(I));
colormap(jet)

% caxis([-0.07 0.1]) % Without Normalization
caxis([-0.012 0.04])
plane_2=plane;
title(['Plane number ', num2str(plane)]);

subplot(1,3,3)
imagesc(squeeze(img(:,:,(2*(plane-229)+50))));
line([0,640],[plane_2,plane_2]);
set(gca,'YDir','normal')
%colormap(gray)

M(i)=getframe(fig);


end


% Write the Video
video=VideoWriter(['C:\Users\robersl\Downloads\video_3.mp4'],'MPEG-4');
video.FrameRate=10;
open(video)
writeVideo(video,M);
close(video)
%%


% for i=1:fileSize
%     char=sprintf('%04d',i);  
%     %h=figure(1);imshow(img(:,:,i),[5000 30000]);pause(0.05)
%     %h=figure(1);imshow(img(:,:,i),[2000 10000]);pause(0.05)
%     imwrite(img(:,:,i),strcat(char,'.png'),'png','bitdepth',16);
% end
