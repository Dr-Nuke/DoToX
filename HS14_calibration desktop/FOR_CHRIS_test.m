
%% setup data and perform FBP reconstruction

close all
clear
clc

load('a');
a=circshift(a,12,2);
trimPixels=40; % number of edge pixels to cut
a=a(trimPixels+1:end-trimPixels,:);
[pixels,projections]=size(a);
stepSize=360/projections;
angles=[0:stepSize:360-stepSize];
tic
disp('Solving with FBP...')
I=iradon(a,angles,'spline','Hann',1,pixels);
toc
I=normalizeMatrix(I);
imshow(I);

%% ART reconstruction

% scale original sinogram as desired
scaling=0.5; % scaling of original sinogram
[rows,cols]=size(a);
a=imresize(a,[round(rows*scaling) round(cols*scaling)]);
% calculate new projection step size
[pixels,projections]=size(a);
stepSize=360/projections;
theta=[0:stepSize:360-stepSize]; % previously calle "angles"

% convert sinogram to column form for use in Ax=b
b=reshape(a,pixels*projections,1);

N=350; % reconstructed image is NxN pixels
p=pixels; % this should be done more properly
d=pixels; % this should be done more properly

% create ART weight matrix
disp('Creating ART weight matrix...')
tic
A=paralleltomoEditRA(...
    N,...
    theta,...
    p,... % p 
    d); % d
    %0); % isDisp
toc

% solve for image using ART
disp('Solving with ART with different K values...')
figure
tic
plotCounter=0;
for K=[5 10 20 50 100 200 400 800 1000]%[2:2:30 35:5:80 90:10:300]
    plotCounter=plotCounter+1;
    [x,info,restart] = AIRtools_SART(A,b,K);
    toc
    I2=reshape(x,N,N);
    I2=normalizeMatrix(I2);
    %I2=adjustScale(I2,2);
    subplot(3,3,plotCounter);
    imshow(I2);
    title(['K = ' num2str(K)])
end

% tic
% disp('Solving with ART with different one K value...')
% K=100;
% [x,info,restart] = AIRtools_SART(A,b,K);
% toc
% I2=reshape(x,N,N);
% I2=normalizeMatrix(I2);
% %I2=adjustScale(I2,2);
% figure;
% imshow(I2);
% title(['K = ' num2str(K)])


% old code  
%imwrite(IART,[folder filename(1:end-4) '_ART_K' num2str(1000+K) '.jpg'],'jpg')
%imshow(IART,[0 1],'Border','Tight');%title('algebraic reconstruction image')
