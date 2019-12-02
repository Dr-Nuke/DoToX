

i=1 % reference
rfname=sprintf('1_corr_%02d.mat',i);
d1=f.loadSingleVariableMATFile(strcat(T.d.DataPath,rfname));


i=11 % flow case, y shifted
rfname=sprintf('1_corr_%02d.mat',i);
d11=f.loadSingleVariableMATFile(strcat(T.d.DataPath,rfname));
%% 
d111=f.fraccircshift(d11,[T.shift(i,1),0,T.shift(i,2)]);
%%

masky=136:348; % spacer area
maskf=30:1460; % all frames
yshiftmask=zeros(size(ysref),'single');
yshiftmask(masky,maskf)=1;
im= squeeze(d111(322,:,:));
imref=squeeze(d1(322,:,:));
[shift,dof]=f.SinoMatchY(im,imref,0,5,21,yshiftmask,-1);
im2=f.fraccircshift(im,[shift,0]);

figure(121)
imshow(imref+yshiftmask,[])

%%
figure(119);clf;
col =[0,1]
ax(1)=subplot(4,1,1);
imshow(squeeze(d1(:,50,:)),col)
title('reference sinogram')

ax(2)=subplot(4,1,2);
imshow(squeeze(d11(:,50,:)),col)
title('flow sinogram')

ax(3)=subplot(4,1,3);
imshow(squeeze(d111(:,50,:)),col)
title('f & x-shifted flow sinogram')

ax(4)=subplot(4,1,4);
imshow(squeeze(d111(:,50,:))-squeeze(d1(:,50,:)),[-0.1,0.1])
linkaxes(ax)
%%




figure(120);clf;
col =[0,1];
bx(1)=subplot(2,2,1);
imshow(imref,col)
title('reference trans sinogram')

bx(2)=subplot(2,2,2);
imshow(im,col)
title('shifted flow trans sinogram')

bx(3)=subplot(2,2,3);
imshow(im-imref,[-0.1,0.1])
title('difference trans flow sinogram')

bx(4)=subplot(2,2,4);
imshow(im2-ref,[-0.1,0.1])
title(sprintf('difference trans flow sinogram x shifted %.3f',shift))
linkaxes(bx)
