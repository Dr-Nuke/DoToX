  clear all;close all
  
  % load image and view it
  I = imread('coins.png');
  figure(1); imagesc(I); colormap(gray); axis equal; axis tight; hold on;
  % (i take imagesc rather than imshow for the axis ticks)
  
  %find the circles & centers & plot
  [centers, radii, metric] = imfindcircles(I,[15 30]);
  viscircles(centers, radii,'EdgeColor','b');
  
  % plotting them works well:
  plot(centers(:,1),centers(:,2));
  
  % now lets modify the image around the first center
  x=round(centers(1,1))-5:round(centers(1,1))+5; %xrange around 1st center
  y=round(centers(1,2))-5:round(centers(1,2))+5; %yrange around 1st center
  
  I(x,y)=0; % put its surrounding to zero
  
  figure(2)
  imagesc(I); colormap(gray); axis equal; axis tight; hold on;
  
  % => fail!
  
  % the plotting has different coordinates:
  plot(x,y,'.');
  
  %%
I = imread('coins.png');
figure()
[a,b]=size(I);
plot([0,a],[0,b]); hold on
imagesc(I);colormap(gray); axis equal; axis tight;
bx=[1,50]
by=[10,100]
I(bx(1):bx(2),by(1):by(2))=0;
imagesc(I');colormap(gray); axis equal; axis tight;
plot([bx(1),bx(2),bx(2),bx(1)],[by(1),by(1),by(2),by(2)])
% => 2 diferent boxes!
