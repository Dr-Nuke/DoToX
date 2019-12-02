
for i=1:4

        rfname=sprintf('2_add_%02d.mat',i);
        d2=f.loadSingleVariableMATFile(strcat(T.d.DataPath,rfname));
        frame(:,:,i)=squeeze(d2(322,:,:));

end

%%

col=[-0.4,0.4];
%ref=squeeze(frame(:,:,1));
ref=ysref;
%b=squeeze(frame(:,:,2));
b=ysim;

% mask=zeros(size(ref));
% mask(136:348,20:1350)=1;
% maskx=136:348;
% masky=20:1350;
mask=yshiftmask;
masky=136:348; % spacer area
maskf=30:1460; % all frames

im=b-ref;
im=im(masky,maskf);

figure(116);clf;
ax(1)=subplot(2,1,1);
imshow(im,col),colorbar;
dif=abs(im);
difff=sum(dif(:));
title(sprintf('original, div %f',difff))

[shift,dof]=f.SinoMatchY(b,ref,0,5,21,mask,-1);
b2=f.fraccircshift(b,[shift,0]);
im2=b2-ref;
im2=im2(masky,maskf);

figure(116);
ax(2)=subplot(2,1,2);
imshow(im2,col),colorbar;
dif=abs(im2);
difff=sum(dif(:));
title(sprintf('shifted by %.3f, div %f',shift,dof))
% guess=[2,T.shift(11,2)]





































