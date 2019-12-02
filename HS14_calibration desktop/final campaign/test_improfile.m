xi=[90,230];
yi=[260,70];
n=100

figure(1);clf;
plane=100;
im=recon(:,:,plane,1);
imshow(im',[0.0001*[-30,40]]);set(gca,'YDir','normal')
hold on
title('kein film')
plot(xi,yi,'r')



[cx1,cy1,c1] = improfile(im,xi,yi,n)

figure(2);clf;
plane=500;
im2=recon2(:,:,plane,6);
imshow(im2',[0.0001*[-30,40]]);set(gca,'YDir','normal')
hold on
title('mit film')
plot(xi,yi,'r')

[cx2,cy2,c2] = improfile(im2,xi,yi,n)

figure(3);clf;
d=sqrt((cx1-cx1(1)).^2+(cy1-cy1(1)).^2);
plot(d,c1,'displayname','ohne film')
hold on

plot(d,c2,'displayname','mit film')
legend()
grid on
xlabel('path')
ylabel('attenuation density')