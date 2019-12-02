%% script to produce bitmaps of the channel cross section
% for robert
close all
clear
% settings
res=0.01; % sampling of the bitmap pixels in mm (i.e. width of a pixel)
domain=30;    % mm size of the image
lfts=0.1:0.1:1.5; % set ofliquid film thickness in mm

%% make an image!
% CAD parameters
r1=2;
r2=3;
r4=4.14;
r5=5.14;
r6=14.604

nx=domain/res+1;  % number of x-pixels
ny=domain/res+1;  % number of y-pixels
xcoord=linspace(-domain/2,domain/2,nx);
ycoord=linspace(-domain/2,domain/2,ny);
[xx,yy]=ndgrid(xcoord,ycoord); % coordinates matrix
P=zeros(nx,ny,length(lfts)); % preallocate the phantom

for k=1:length(lfts)
    disp(sprintf('getting lft %2.1f',lfts(k)))
    % generate phantom, focus on 1 quarter
    Ph = f4_ChannelPhan2_w_film(lfts(k));
    Ph = f42_RotPhan(Ph,pi/4)
   
    % air is where...:
    air1=xx.^2+yy.^2>(r6)^2; % outer circle
    air2=sqrt((xx-Ph.c(5,1)).^2+(yy-Ph.c(5,2)).^2)<r1; % small cirlce outside
    
    s=4; %segment number
    air31=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
    air32=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s));
    air33=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r2;
    air3=air31 & air32 & air33;
    
    s=6;
    air4=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
    air4(atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s)))=0;
    air4(sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r2)=0;
    % air4=air41 & air42 & air43;
    
    air5=(yy)>Ph.b(18);
    air5(xx<Ph.Lx2(18,1))=0;
    %air5(yy<Ph.Lx1(18,2))=0;
    air5(xx>Ph.Lx1(18,1))=0;
    
    air=air1|air2|air3|air4|air5;
        
    s=16;
    d2o=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r4;
    d2o((yy)>Ph.b(17))=0;
    s=2;
    d2o((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)-2*pi) &...
        atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s)-2*pi)) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r1)=0;

    s=1;
    liq=zeros(nx,ny);
    
    s=13;
    liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
        atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lfts(k)))=1;
    
    s=1;
    liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
        atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lfts(k)))=1;
    
    s=15;
    liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)-2*pi) |...
        atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s))) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>Ph.r(s) &...
        sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<(Ph.r(s)+lfts(k)))=1;
    
    %assemble mask
    mask=air+2*d2o+4*liq;
    
    % mirror that 1/8th to make a quarter
    for i=ceil(nx/2):nx
        for j=ceil(nx/2):i
            mask(i,j)=mask(j,i);
            
        end
    end
    % mirror quadrants
    mask(1:ceil(nx/2),ceil(nx/2):end)=...
        flipud(mask(ceil(nx/2):end,ceil(nx/2):end));
    % make full mask
    mask(:,1:ceil(nx/2))=...
        fliplr(mask(:,ceil(nx/2):end));
    
    % get vapor area and set it to "3"
    [tempmask,n]=bwlabel(~im2bw(mask));
    mask(tempmask==tempmask(ceil(nx/2),ceil(nx/2)))=3;
    
    % archive the mask
    P(:,:,k)=mask;
    figure(123);clf;
    imshow(P(:,:,k)',[]);set(gca,'Ydir','normal')
end

save('P','P')















