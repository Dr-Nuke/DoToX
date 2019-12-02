figure(23)
clf


x=[350,2115];
y=[0.05,1]
xq=[500,1000,1500,2000];
mm=interp1(x,y,xq);

cmax=0.0096;
for i=1:4
%imshow(squeeze(recontot(:,:,k))',[],'Colormap',hsv); axis equal; axis tight;set(gca,'YDir','normal');

subplot(2,2,i)
ax1=imshow(squeeze(recontot(:,:,xq(i)))',[0,cmax]); axis equal; axis tight;set(gca,'YDir','normal'); 
title(sprintf('%4.2f mm LFT',mm(i)))
end

% subplot(2,2,2)
% ax1=imshow(squeeze(recontot(:,:,1000))',[0,cmax]); axis equal; axis tight;set(gca,'YDir','normal'); 
% 
% subplot(2,2,3)
% ax1=imshow(squeeze(recontot(:,:,1500))',[0,cmax]); axis equal; axis tight;set(gca,'YDir','normal');
% 
% subplot(2,2,4)
% ax1=imshow(squeeze(recontot(:,:,2000))',[0,cmax]); axis equal; axis tight;set(gca,'YDir','normal');