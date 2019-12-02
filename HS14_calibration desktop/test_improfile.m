close all
%im=f_KernelGen(19,18,14);

im=f_normalize(im2double(fitsread('CL3_0040.fits')'));
x=[244,362];
y=[246,361];

x=[200,600];
y=[300,700];
len=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2)


imbo3(im,1)
hold on
plot(y,x,'Color','r','LineWidth',2)

figure(2)
[cx,cy,c,xi,yi]=improfile(im,x,y,15*len,'bilinear');
plot(c,'.')
length(c)








