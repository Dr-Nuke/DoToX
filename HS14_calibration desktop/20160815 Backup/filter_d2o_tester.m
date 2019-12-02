close all
path='C:\data\tomo_HS14\02_rawdata\D2O\D2O_0012.fits';
t=5000
im=im2double(fitsread(path))';
im2 = medfilt2(im,[3,3]);
im3=im-im2;

immin=min(im(:))
immax=max(im(:))


n=5; %number of subplots

mymap = [0 0 0; 0.5 0.5 0.5; 1 1 1]; %colormap

figure(1)
subplot(1,n,1)
ax1=imagesc(im');colormap(gray); axis tight;set(gca,'YDir','normal');
caxis([0 60000]);

subplot(1,n,2)
ax2=imagesc(im2');colormap(gray); axis tight;set(gca,'YDir','normal');
caxis([0 60000]);

subplot(1,n,3)
ax3=imagesc(im3');colormap(gray); axis tight;set(gca,'YDir','normal');

caxis([-t t]);

im4=im.*(abs(im2-im)>0.05*im);

imf=im;imf(abs(im2-im)>0.05*im);


% 
% subplot(1,n,4)
% ax4=imagesc(im3');colormap(gray); axis tight;set(gca,'YDir','normal');
% caxis([-t t]);
% 
% figure(2)
% ax3=imagesc(im3');colormap(gray); axis tight;set(gca,'YDir','normal');
% caxis([-t t]);
% colormap(mymap)