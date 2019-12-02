close all
i=1500;
im=f_loadimage(i);

tic;
%small spots
thresh=0.02;
im2=f_MediFilter1(im,thresh);
t(1)=toc

%camera lines
im3=f_CamLineFilt(im2);
t(2)=toc


% custom median
x1=600;x2=1200;y1=300;y2=1000;
kernel=f_KernelGen(9,7.5);
mode='median';
tmp=f_MediCustom(im3,kernel,mode);

%%
t(3)=toc
im4=im3;
% make the difference
diff=im3-tmp;

% correct output
threshold=0.06
im4(abs(diff)>threshold*im4)=tmp(abs(diff)>threshold*im4);



 close all
% figure('Name','raw','units','normalized','outerposition',[0 0 1 1])
% imbo4(im');axis tight;axis equal;
% %axis([x1,x2,y1,y2])
% 
% figure('Name','smallspots','units','normalized','outerposition',[0 0 1 1])
% imbo4(im2');axis tight;axis equal;
% %axis([x1,x2,y1,y2])

figure('Name','lines','units','normalized','outerposition',[0 0 1 1])
imbo4(im3');axis tight;axis equal;

figure('Name','custom median','units','normalized','outerposition',[0 0 1 1])
imbo4(im4');axis tight;axis equal;

figure('Name','custom diff','units','normalized','outerposition',[0 0 1 1])
imbo4((im4.*(diff>threshold*im4))');axis tight;axis equal;

