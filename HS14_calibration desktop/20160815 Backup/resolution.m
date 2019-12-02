figure(4)
cla
im_res=fitsread('reso.fits');
im_res=imrotate(im_res,1.57,'nearest','crop');
figure(1)
imbo4(im_res)
hold on
x=106;
y=161;
h=2210;
w=1991;
%plot([,x],[1,2560],'r')
rectangle('Position',[x,y,w,h],'EdgeColor','r')

res_ud=h/100 % pixel/mm
res_lr=w/90  % pixel/mm