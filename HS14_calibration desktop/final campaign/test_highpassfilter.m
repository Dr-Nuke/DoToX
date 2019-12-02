% make data
im=mean(d,3);
h=figure(500);clf
set(h,'position',[-275 0 1680+275 1050])

spacermask=ones(size(im));
spacermask(50:270,150:330)=0;


% plot data
crange=[0.5 0.65];
thresh=0.006;
nplots=7;
height=1;
width=1/nplots;
ax(1)=subplot('position', [0 0 width height]);
imagesc(im');set(gca,'YDir','normal')
colormap(gca,'hsv')
caxis(crange)
axis equal
axis off
title('original')

% median filter
med=medfilt2(im,[3 3]);
mask1=(abs(im-med)>thresh);

ax(2)=subplot('position', [width 0 width height]);
imagesc(mask1');set(gca,'YDir','normal')
colormap(gca,'gray')
axis equal
axis off
title('median mask')


imcorr1=im;
imcorr1(mask1)=med(mask1);

ax(3)=subplot('position', [2*width 0 width height]);
img1=imagesc(imcorr1');set(gca,'YDir','normal')
colormap(gca,'hsv')
caxis(crange)
axis equal
axis off
title('median applied')

%gauss filter
 [x, y]=meshgrid(round(-1):round(1), round(-1):round(1));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 f=f./sum(f(:));
gaussim=conv2(imcorr1,f,'same');


ax(4)=subplot('position', [3*width 0 width height]);
imagesc(gaussim');set(gca,'YDir','normal')
colormap(gca,'hsv')
axis equal
axis off
title('median mask')

% ax(5)=subplot('position', [3*width 0 width height]);
% imagesc(imcorr2');set(gca,'YDir','normal')
% colormap(gca,'hsv')
% caxis(crange)
% axis equal
% axis off
% title('gauss applied')

% camera spot removal

k=mt.KernelGen(11,3.5,0.6);



%[x,y]=ndgrid(1:size(im,1),1:size(im,2))




linkaxes(ax)
%%

h2=figure(501);clf
set(h2,'position',[0 0 1680 1050])
% plot data
crange=[0.5 0.65];
thresh=0.01;
nplots=6;
height=1;
width=1/nplots;
bx(1)=subplot('position', [0 0 width height]);
imagesc(im');set(gca,'YDir','normal')
colormap(gca,'hsv')
caxis(crange)
axis equal
axis off
title('original')

rs=[1 2 3 4 5]
threshs=0.0002*logspace(1,2,6);
for i=2:6
   
   bx(i)=subplot('position', [(i-1)*width 0 width height]); 
   %k=mt.KernelGen(11,rs(i-1),0.6);
   imcorr3=conv2(im,k,'same');
   imagesc((im'-imcorr3')>threshs(i));set(gca,'YDir','normal')
   colormap(gca,'gray')
   %caxis([-0.01,0.01])
   axis equal
   axis off
   title('ring gauss applied')
       

end

linkaxes([ax bx])

%%
badim=im(120:130,470:480);
badimsize=size(badim);
[x,y]=meshgrid(1:badimsize(1),1:badimsize(2));

mask=~mt.KernelGenBin(11,5,0);
badimring=badim(mask);
xring=x(mask);
yring=y(mask);






