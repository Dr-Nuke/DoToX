%% create a phantom


addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

%phantom image pixel size. must be odd
phan_size=201;

%check for odd-ness
if mod(phan_size,2)==0
    error('the phantom size is not odd')
end

%initialize phantom
phan=zeros(phan_size);

%create coordinates 
[x,y]=ndgrid(1:phan_size,1:phan_size);
%x=x-(ceil(phan_size/2));
%y=y-(ceil(phan_size/2));

%mask
xx=[floor(phan_size/4),floor(3*phan_size/4)];
yy=floor(phan_size/4):floor(3*phan_size/4);

%two parallel lines, sinus intensities
phan(xx,yy)=sin(10*x(xx,yy)+6*pi*y(xx,yy)/max(max(y(xx,yy))));

%two circles
%phan((x.^2<(10^2-y.^2,:)=1



% figure(1)
% clf
% imbo4(phan)


