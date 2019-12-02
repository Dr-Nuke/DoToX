%script to find the thickness of the liquid film (bright areas) on the 
% coolant rods (dark/ not seen)


% file names and paths
clc
clear all
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

block=f3_loadSingleVariableMATFile('C:\data\Tomo_HS14\processed\block_7rec_cl3_.mat');
%%
plane=300;
cl33=imresize(squeeze(block(:,:,plane)),5);
%cl33=imrotate(cl33,-5,'bilinear');
cl333=f_normalize(cl33);

%create black/white-image (helps the algorhythm)
clbw=im2bw(cl333,0.6); %empirically chosen threshold of 0.6 for greylevel to bw coonversion
%clbw=imcomplement(clbw); %invert the bw image

res=22.1;       % pixel/mm, camera characteristic
r_rod=5.14;     % physical radius of the coolant rods
dr=3;           % deviation around the radius in pixel
r_min=round(r_rod*res-0.5*dr); %5.14 mm physical radius of the coolant rods
r_min=round(112);      % empirical test shows the above calculation underestimates
S=0.99;         % Sensitivity threshold
o =40;          % overlap


%% find centers & radii
[c,rad,m]=f_hugh_4(clbw,r_min,dr,S,o,0.6); %have this file in the same folder
c2=fliplr(c);
disp([c2,rad,m])

%% plot result
figure(1);
clf
imshow(cl33',[]);colormap(gray);axis equal; axis tight;set(gca,'YDir','normal');
hold on
scatter(c2(:,1),c2(:,2),'gx') %centers
viscircles(c2, rad,'EdgeColor','g','LineStyle',':','LineWidth',1,'EnhanceVisibility',0);
axis([0 749 0 749]); axis equal




%% find angles and record improfiles
figure(1)

hold on
da=135; % angel incrment
n=da+1; % number of angles
r2=150; % length of path in pixel (f_LinGen) and number of samples (improfile)

for i=1:4,
    j=mod(i,4)+1;
    disp([i,j]);
    a_start=radtodeg(atan2(c2(j,2)-c2(i,2),c2(j,1)-c2(i,1)))-(da-90)/2;
    [px,py] = f_LinGen(c2(i,:), a_start,da,n,r2);
    for j=1:n,
        plot([c2(i,1),px(j)],[c2(i,2),py(j)]) %plot beams
        [cx,cy,c,xi,yi]=improfile(cl33',[c2(i,1),px(j)],[c2(i,2),py(j)],r2); % generate improfiles
        c_tot(i,j,:)=c;
        cx_tot(i,j,:)=cx;
        cy_tot(i,j,:)=cy;
    
    end
    
    pxx(i,:)=px;
    pyy(i,:)=py;
end
plot([c2(:,1);c2(1,1)],[c2(:,2);c2(1,2)],'r')
%axis([0 749 0 749]); axis equal

%% plot line profiles
figure(2)
clf
imagesc(cl33');
set(gca,'YDir','normal'); 
colormap('gray')
xlabel('x');
ylabel('y');
zlabel('z');
col=hsv(n)
 
for i=1:4



    hold on
    for j=1:n
        plot3(squeeze(cx_tot(i,j,:)),...
              squeeze(cy_tot(i,j,:)),...
              squeeze(c_tot(i,j,:)),'color',col(j,:))
    end
        


end


 
%% do the line profiles





