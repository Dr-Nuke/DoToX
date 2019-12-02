clear all
clc
format compact
close all

addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
% load film movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1134h/imc2.mat')

% load empty movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1341h/imc1.mat')
%%
imc2=imc2(:,:,1:1455);
yline=45;

a=double(squeeze(imc1(:,yline,:)));
b=double(squeeze(imc2(:,yline,:)));


figure(1);clf;
subplot(1,3,1)
imshow(a',[]);set(gca,'YDir','normal')
title('film sinogram')
subplot(1,3,2)
imshow(b',[]);set(gca,'YDir','normal')
title('empty sinogram')

subplot(1,3,3)
diff=a'-b';
imshow(diff,[]);set(gca,'YDir','normal')

d=sum(abs(diff(:)))
title(sprintf('difference sinogram: %.0f deviation',d ))

%% check xcov, frame dimension

for i=2:size(a,1)-1
    [c,lags]=xcov(a(i,:),b(i,:),'coeff');
    [~,ind]=max(c);
xshift(i)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
end

% cut out outliers
meanfshift=mean(xshift(abs((xshift-mean(xshift))/mean(xshift))<0.2))

figure(2);clf;
plot(xshift)
ylim([-12.6,-12.15])
xlabel('x position [pixel]')
ylabel('lag')
grid on
title(sprintf('cross correlation based lag. %.2f mean',meanfshift))


%% shifting the image, Frame direction
cbox=[-0.1,0.1];
%b2=f4_ShiftY(b,meanxshift);
b2=f4_ShiftXY(b,[0,meanfshift]);
figure(3);clf;
ax(1)=subplot(1,3,1);

imshow(diff,cbox);set(gca,'YDir','normal')
title(sprintf('original difference (%.0f dev.)',d ))

diff2=a'-b2';
ax(2)=subplot(1,3,2)
imshow(diff2,cbox);set(gca,'YDir','normal')
d2=sum(abs(diff2(:)));
title(sprintf('diff. after Frame-shift (%.0f dev.)',d2 ))

%% x shift
for i=1:size(a,2)
    [c,lags]=xcov(a(:,i),b2(:,i),'coeff');
    [~,ind]=max(c);
yshift(i)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
end

meanxshift=mean(yshift(abs((yshift-mean(yshift))/mean(yshift))<0.2));

figure(4);clf;
plot(yshift,1:length(yshift))
ylabel('Frame number')
xlabel('lag')
grid on
title(sprintf('cross correlation based x-lag. %.2f mean',meanxshift))

b3=f4_ShiftY(b,meanxshift);
figure(3)


%% do the x shift
b3=f4_ShiftXY(b,[meanxshift,meanfshift]);
%b3=f4_ShiftXY(b,[0,0]);
figure(3);
subplot(1,3,3)



diff3=a'-b3';
ax(3)=subplot(1,3,3)
imshow(diff3,cbox);set(gca,'YDir','normal')
d2=sum(abs(diff3(:)));
title(sprintf('difference after both shifts (%.0f dev.)',d2 ))

%% iterate 2d shift grid
figure(5);clf
sxb=linspace(-20,20,20); %basis iteration grd
sfb=linspace(-20,20,20);
sx=sxb+meanxshift;
sf=sfb+meanfshift;

double(squeeze(imc2(:,400,:)));
for h=1:3

diff2d=zeros(length(sx),length(sf));

for i=1:length(sx) % x dimension
    for j=1:length(sf) % frame dimension
        k=j+(i-1)*length(sx);
        f_BoCount(k,20,10,5)
        im=f4_ShiftXY(b,[sx(i),sf(j)]);
        diff2d(i,j)=sum(abs(im(:)-a(:)));
    end
end

% optimize for region
[M,I] = min(diff2d(:));
[I_row, I_col] = ind2sub(size(diff2d),I);
subplot(1,3,h)
imagesc(diff2d','XData', sx, 'YData', sf);set(gca,'YDir','normal');colorbar
title(sprintf('2D deviation by pixel shift\n iteration %d, min dev. %.0f',...
    h,diff2d(I_row,I_col)))

sxb=sxb/10;
sfb=sfb/10;
sx=sxb+sx(I_row);
sf=sfb+sf(I_col);
end

        
        
%% hard code check
im=f4_ShiftXY(b,[sx(I_row),sf(I_col)]);
diff4=a'-im';
d3=sum(abs(a(:)-im(:)));

match=[sx(I_row),sf(I_col)];
save('match','match');
       
        
figure(3);
subplot(1,3,3)
imshow(diff4,cbox);set(gca,'YDir','normal')
title(sprintf('difference after both shifts (%.0f dev.)',d3 ))
        
        
%% check if this is true for other heights

sxb=linspace(-5,5,10); %basis iteration grd
sfb=linspace(-5,5,10);

planediff=zeros(1,size(imc2,2));
planeshift=zeros(size(imc2,2),2);

parfor g=1:size(imc1,2)
disp(g)
b=double(squeeze(imc2(:,g,:)));

    for h=1:2
    sx=sxb/(10^(h-1))+match(1);
    sf=sfb/(10^(h-1))+match(2);    

    diff2d=zeros(length(sx),length(sf));

    for i=1:length(sx) % x dimension
        for j=1:length(sf) % frame dimension

            im=f4_ShiftXY(b,[sx(i),sf(j)]);
            diff2d(i,j)=sum(abs(im(:)-a(:)));
        end
    end

    % optimize for region
    [M,I] = min(diff2d(:));
    [I_row, I_col] = ind2sub(size(diff2d),I);

    end % iterations

planediff(g)=M;
planeshift(g,:)=[sx(I_row),sf(I_col)];


end

%%




























