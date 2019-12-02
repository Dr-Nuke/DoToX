
% finds the pixel resolution
%% en Files
clear all
close all
clc


filepath=('C:/data/20161216 Campaign 5 Robert Z/11h34 Tmax 110kV 5mA.seq');
fileInfo = dir(filepath);
imsize=[1024 640];
header_size=2048;
fileSize = fileInfo.bytes;
fileSize=floor((fileSize-header_size)/imsize(1)/imsize(2)/2);
fid=fopen(filepath);
hed = fread(fid,header_size);%the header size migth need to be adjusted depending on image settings
%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
img=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16')));
fclose(fid);
img=rot90(reshape(img,[1024 640 fileSize]));

%% find all the relevant frames and dispose the rest
diff=img(:,:,1:end-1)-img(:,:,2:end);
diffind=squeeze(max(max(diff)));
figure(1);clf;
plot(diffind,'x-')
xlabel('frame #');ylabel('abs(mean(diff))');title('frame difference plot')

%manual check
startframe=27;
endframe=1422;

figure(2);clf
imax=6;
for i= 1:imax
    subplot(imax,2,2*i-1)
    imshow(squeeze(diff(:,1:100,startframe-2+i))',[0 10000]);set(gca,'YDir','normal')
    title(sprintf('frame %d - frame %d',startframe-2+i,startframe-2+i+1 ))
    
    subplot(imax,2,2*i)
    imshow(squeeze(diff(:,1:100,endframe-4+i))',[0 10000]);set(gca,'YDir','normal')
    title(sprintf('frame %d - frame %d',endframe-4+i,endframe-4+i+1 ))
end

% check again manually
startframe=28;
endframe=1423;
clear diff ;
imgcut=img(:,:,startframe:endframe);
clear img

%% find the dose and beam correction area
figure(3);clf
imbo4(imgcut(:,:,1))
hold on

dose=[[5,427,429];... % x start value
    [162,630,627];...      % x end value
    [100,100,700];...      % y start value
    [1010,229,1014]];         % y end value

doseimg=nan(size(imgcut(:,:,1)));
for i=1:3
    doseimg(dose(1,i):dose(2,i),dose(3,i):dose(4,i))=...
        imgcut(dose(1,i):dose(2,i),dose(3,i):dose(4,i),1);
end


%check if this makes sense
figure(4);clf;
imbo4(doseimg)
title('this shows the dose correction area bases')
xlabel('x')
ylabel('y')
%check if nothings in these areas
figure(5);clf;
dosemask=uint16(~isnan(doseimg));
dosemask3=uint16(repmat(dosemask,[1,1,size(imgcut,3)]));
dosemask3=dosemask3.*imgcut;

for i=1:size(imgcut,3)
   a=dosemask3(:,:,i);
   dosemin(i)=min(a(a>0));
end
plot(dosemin)
ylabel('minimum image value of masked frames')
xlabel('frame number')
title('if this plot has dips, then some object rotate into the dose mask')
clear dosemask3

%% do the beam nonuniformity correction (thanks to lukas)

% Create point (x,y,z) of the background data for the 2D-fit
meshsize=size(dosemask);

[x,y]=ndgrid(1:meshsize(1),1:meshsize(2));
x(dosemask==0)=[];
y(dosemask==0)=[];


% find the fits
for i=1:size(imgcut,3)
    f_BoCount(i,50,10,5)
    im=squeeze(imgcut(:,:,i));
    im(dosemask==0)=[];
    f{i}= fit( [x',y'],im', 'poly22' );
    
end
    
%% check if the fits are good
figure(5);clf;
for i = round(size(imgcut,3)*rand(1,10));

im=squeeze(imgcut(:,:,i));
im(dosemask==0)=[];
rf =10; % reduction factor

plot(f{i},[x(1:rf:end)',y(1:rf:end)'],im(1:rf:end)')
hold on
%plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
title(sprintf('nun-uniformity data and according fit of frame %d\n together with the fit',i))
xlabel('x')
ylabel('y')
zlabel('counts')
pause(1)
end

%% so now correct for it

[x,y]=ndgrid(1:meshsize(1),1:meshsize(2));
unibeam=zeros(size(imgcut),'single');

for i=1:size(imgcut,3)
    f_BoCount(i,50,10,5)
    unibeam(:,:,i)=feval(f{i},x,y);
    
end

% apply the correction

savefast('imgcut','imgcut');
imgcut=single(imgcut)./unibeam;
clear unibeam

%%






    %%
figure(2)						% Check if data looks good
plot(z)
f = fit( [x,y],z, 'poly22' );	% Do the fit
figure(3)						% Show the fit
plot(f,[x(1:100:end),y(1:100:end)],z(1:100:end))

% Create Reference image from fit
for i=1:size(img_ref,1)
    for j=1:size(img_ref,2)
        img_ref(i,j)=f(i,j);
    end
end
%%
% Normalize each frame with the Background
for i=1:size(img,3)
    img(:,:,i) = img(:,:,i)./img_ref;
end





%% find area for resolution check
testimg=squeeze(img(:,:,1));
imshow(testimg',[]);set(gca,'YDir','normal')
title('sample image for xray camera resolution');
hold on
rpp=[[184,184];[220,630]]; % resolution profile points

x=443;
y=224;
h=420;
w=22;
rp=[x y w h];

rpp=[x,x,x+w,x+w x; y,y+h, y+h, y,y];

rectangle('Position',rp);
figure(68);clf;
subplot(1,1,1);hold on;
col=hsv(w);


for i=1:w
    xi=x*ones(1,h+1)+i;
    yi=[y:y+h];
    [cx,cy,c] = improfile(testimg',xi,yi,h+1);
    plot3(cx,cy,c,'color',col(i,:));
end

xlabel('x pixels');
ylabel('y counts');
title('counts');
%%


%% make a mean 
figure()
imshow(testimg(x:x+w,y:y+w),[]);set(gca,'YDir','normal')
resprofile=mean(testimg(x:x+w,y:y+h),1)
figure(67)
plot(y:y+h,resprofile,'x-')
title('mean resolution profile')

% manual write down of rising & falling edge midde point

fe=[238,322,406,489,573];
re=[280,364,447,531,611];
dif=re-fe;
res=dif/5

% rsult is 8.2 pixel per mm





