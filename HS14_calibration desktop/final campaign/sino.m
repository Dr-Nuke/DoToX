clear all
close all

%% test of the centering
s=301; % rough width of my sinograms

CREATE_VIDEO = false;

% % phantom:
% phan=zeros(s);
% px=120; % impulse position
% phan(px,px)=1;
% phan(200,200)=1;
 ang=360;  % angular range
 nang=1357; % number of angles

% phantom is kindergarden, we use data!
load('filmsino')

sinogram=filmsino(160:490,:);
k180=ceil(size(sinogram,2)/2);

% FIND THE CENTERINGSHIFT sry caps

[c,lags]=xcov(sinogram(:,1),flipud(sinogram(:,k180)),'coeff');
 [~,ind]=max(c);
centershift=(lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
    (log(c(ind-1))-2*log(c(ind))+log(c(ind+1))))/2;
sinogram=fraccircshift(sinogram,-centershift);

ss=size(sinogram);

sinogram(98,:)=1;
% angles list
angles=linspace(1,ang,nang+1);
angles(end)=[];

% % show phantom
% fig=figure(8);clf;
% fig.Position(1:2)=[0,0];
% imagesc(phan');set(gca,'YDir','normal')colorbar
% setfig(fig)
% title('phantom')

% % do the iradon
% [sinogram,xp]=radon(phan,angles);

% % crop sinogram (dirty joke haha
% cs=(size(sinogram,1) - s)/2; %CropSize
% sinogram=sinogram(cs+1:end-cs,:);


% show sinogram
fig=figure(11);clf;
fig.Position(1:2)=[500,500];
imagesc(sinogram');set(gca,'YDir','normal')
title('sinogram')
setfig(fig)

%initialize the result
recon=zeros(ss(1));

% s2=length(xp);
h=figure(13241);clf;

set(h,'Position',[-500         10        2100        523])
sub(1)=subplot(1,3,1);

sub(2)=subplot(1,3,2);

sub(3)=subplot(1,3,3);

title('reconstruction');

sinogram_gpu = gpuArray(sinogram);
sinogram_plot = ones(size(sinogram));

if CREATE_VIDEO
    myvideo = VideoWriter('bo_movie.avi');
    open(myvideo);
end

for ang=1:length(angles)
   % creat backprojection
   bp=repmat(sinogram(:,ang),[1,ss(1)]); % one projection smeared out
   bp=imrotate(bp,angles(ang),'bilinear','crop');
   recon=recon+bp/nang;
   
   if mod(ang,5) == 0
       % subplot 1
       sinogram_plot(:,1:ang) = sinogram(:,1:ang);
       imagesc(sub(1),sinogram_plot')
       set(sub(1),'Ydir','normal');
       
       % subplot 2
       rec2=recon/max(recon(:));
       
       imagesc(sub(2),rec2');set(sub(2),'Ydir','normal');
       colormap(gray)
       caxis([-0.005,0.01])
       
       % subplot 3
       recon2 = iradon(sinogram_gpu(:,1:ang),angles(1:ang),'Hann',1,ss(1));
       set(gca,'Ydir','normal');
       subplot(1,3,3)
       imshow(recon2',[-0.005,0.01])
       
%        colorbar;
%        pause(0.1)
        drawnow;
        if CREATE_VIDEO
            frame = getframe(h);
            writeVideo(myvideo, frame);
        end
   end
end

if CREATE_VIDEO
    close(myvideo);
end


figure(16586);clf
imshow(iradon(sinogram(:,1:ang),1:ang,'spline','Hann',1,ss(1))',[])
set(gca,'Ydir','normal')



 
function setfig(fid)
    xlabel('b')
    ylabel('v')
    colormap(gray)
    colorbar
end




