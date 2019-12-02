% this script tests the circle finding
% load('recon')

i=3;
k=400;
imbinthresh=[0.0003,0.00025,0.00022]
im=squeeze(recon(i,:,:,k));
imbw=im2bw(f_normalize(im),0.55); 
imbw2=~imbinarize(im,imbinthresh(i));

[c,rad,m]=f_hugh_4_2(imbw2,F.r_min,F.dr,F.S,F.o,0.4); %have this file in the same folder
            



figure(122);clf
imshow(im',[]);set(gca,'YDir','normal')

figure(123);clf
imshow(imbw2',[]);set(gca,'YDir','normal')

figure(124);clf
histogram(im,100);
