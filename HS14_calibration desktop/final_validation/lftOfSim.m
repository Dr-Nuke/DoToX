%% post process the micha-matches
% ratio refers to the ratio between old and new tomos, i.e. 1.4700898e-04
% robertfactor is the phantom-to-tomo factor of film thickness
fixfac=T.Rec.detPitch/ratio/robertfactor;

% F the film struct
load('T')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/apriori/')
G.pp.loadpath='T:\cbolesch\res1'; % the path for loading the recons
G.pp.casstr=8;
G.pp.plastr=11:14;
[G.pp.flist,G.pp.cases,G.pp.planes,G.raw]=G.GenRecFileList(F,200);
G.raw.maxp=max(G.raw.nplane)


%%
G.size=size(rec);
G.Nfilt=G.size(3);%G.size(1);
% F means Film, varia
G.r_rod=5.14;
G.res=8.2; %manual fit
G.dr=5;                   % deviation around the radius in pixel
%G.r_min=round(G.r_rod*res-0.5*G.dr);  %physical radius in pixel, with tolerance
G.r_min=40;             % empirical check shows 
G.S=0.99;                  %sensitivity threshold
G.o =10;                  % overlap
G.dh=0.8473;         % hydraulic diameter in cm

% seems unused: G.n_path=round((90/360)*2*pi*G.r_rod*res);  % number of sample paths per pin
G.r_path=60; % path length in pixel +manual fix

G.n_pins = 4; % number of fuel pins
G.d_angle=180; % angel range
G.realangle=174.83; % manual check in inventor
G.n_angles=G.d_angle+1; % number of angles
G.h=G.size(3); % stack height

G.phanLFT=lfts;
G.att=xs; %function xs for attenuation coefficients

% find centers

%preallocate
%G.c=zeros(G.Nfilt,G.size(3),4,2);  % center [x,y] coordinates of circles
%G.rad=zeros(G.Nfilt,G.size(3),4);  % radius values of the circles
%G.m=zeros(G.Nfilt,G.size(3),4);    % confidencey of this circle fit
G.a0=zeros(G.Nfilt,G.size(3),4); % first angle for each circle
%G.ProfEndPnt=zeros(G.Nfilt,G.size(3),G.n_pins,G.n_angles,2); % end points of the integration paths
%G.cfit=zeros(G.Nfilt,G.size(3),4,2); % fittet centers
%G.ProfAngles=zeros(G.Nfilt,G.size(3),4,G.n_angles);
%G.Profs=zeros(G.Nfilt,G.size(3),4,G.n_angles,G.r_path+1);
%G.bgHeat=zeros(G.Nfilt,G.h,G.n_pins);
G.profs=zeros(Nfilt,G.h,G.n_pins,G.n_angles,G.r_path+1)
%
%G.centeringrange=1:1024; %range used for fitting the center locations
%G.ImBinThresh=[0,0.0003,0.00025,0.0003]; % individual thesholds for the ...
% bw conversion per case


% find rod centers via CAD
%%
fh=figure(312);clf
imshow( rec(:,:,1,2)',[])
% set(fh,'units','normalized',...
%     'outerposition',[0 0 1 1]);
set(gca,'YDir','normal')
hold on
colorbar
axis on
%
G.cen.prot=45; % deg
G.cen.pzoom=7.9; % factor
G.cen.px=135;
G.cen.py=135;
G.cen.pcenid=[61 13 29 45]; % counter-clockwise, starting east


G.cen.P=f42_ChannelPhanTransform();
G.cen.P2=f42_RotPhan(G.cen.P,deg2rad(G.cen.prot));
G.cen.P2=f42_PhanZoom(G.cen.P2,G.cen.pzoom);
G.cen.P2=f42_PhanMove(G.cen.P2,G.cen.px,G.cen.py);
f42_PlotPhan(G.cen.P2,gca,[1 0 0],[1 4]);

G.cen.cen=G.cen.P2.c(G.cen.pcenid,:);
%%
for i =G.cen.pcenid;%1:length(P2.c)
    plot(G.cen.P2.c(i,1),G.cen.P2.c(i,2),'ob');
    %text(G.cen.P2.c(i,1),G.cen.P2.c(i,2),num2str(G.cen.P2.r(i)),'FontSize',20);
end
%%
% find paths
[G.Pa.ProfEndPnt,G.Pa.ProfAngles]=lft.FindProfilePathsM(G);

%% check path enpoints
for i=1:4
    plot(squeeze(G.Pa.ProfEndPnt(i,:,1)),squeeze(G.Pa.ProfEndPnt(i,:,2)),'.m')
end
save('F','F')
%%

% obtain the profiles
for i = 1:G.raw.ncas
    G.Profs(i,:,:,:,:)=lft.FindProfilesM(G,squeeze(rec(:,:,i,:)));
end

%%
%check a few profiles

figure(1347);clf;
cas=2;
pla=10;
pin=1;
%imshow(rec(:,:,cas,pla),[]);

[xx,yy]=ndgrid(1:269,1:269);
    zz=zeros(269,269);
su=surf(xx,yy,zz,rec(:,:,cas,pla),'edgecolor','none')
colormap(gray)
angs=1:10:G.n_angles;
nang=length(angs);
col=hsv(nang*4);
hold on
for pin=1:4
for ang2=1:nang;
    ang=angs(ang2);
    x=linspace(G.cen.cen(pin,1),G.Pa.ProfEndPnt(pin,ang,1),size(G.profs,5));
    y=linspace(G.cen.cen(pin,2),G.Pa.ProfEndPnt(pin,ang,2),size(G.profs,5));
    plot3(x,y,squeeze(G.Profs(cas,pla,pin,ang,:)),'color',col(ang2+(pin-1)*nang,:))

end
end
%%
figure(2356);clf;
angs=1:20:G.n_angles;
nang=length(angs);
col=hsv(nang);
cas=2;
tstr={'large flow','small flow'};
for sb=1:2
    ax=subplot(1,2,sb);
    hold on
    for ang2=1:nang;
        ang=angs(ang2);
        plot(squeeze(G.Profs(sb+1,pla,pin,ang,:)-G.Profs(1,pla,pin,ang,:))...
            ,'color',col(ang2,:),'displayname',sprintf('%d°',ang));
    end
    if sb==1
        ylabel('signal')
    end
    ylim(1e-5*[-0.5,2.5])
    xlabel('path length')
    title(tstr{sb})
    grid on
end
legend('location','best')




%%
% 1st approach: sum up the integrals
G.LFT=lft.LFT(F);
save('F','F')
%%


for cas=[2,3,4];
fh=figure(5683+cas);clf;
% set(fh,'units','normalized',...
%     'outerposition',[0 0 1 1]);



for sb=1:4
    subplot(1,4,sb);
    im=convn(squeeze(G.LFT(cas,:,sb,:)-G.LFT(1,:,sb,:))*fixfac,mt.ElliKernelGen(1,5,1,2),'same');
    imagesc(im)
    caxis(1e-5*[0,18]*fixfac)
    title(sprintf('[%d %d]',cas,sb))
    
    set(gca,'YDir','normal')
    set(gca,'XDir','reverse')
    %xticks([1 91 81]);
    yticks([]);
    xticks([1 90 181])
    xticklabels({'90°','0°','-90°'})
    
end
colorbar
end

%% profiles
cas=2;
pin=2;

planes=[2 6:5:max(G.raw.nplane)]
offset=1e-4;
angu=linspace(-G.d_angle/2,G.d_angle/2,G.n_angles);

figure(4572);clf
ax=gca;
hold on


for i=1:length(planes)
ph(1)=plot(angu,flipud(squeeze(G.LFT(2,planes(i),pin,:)-G.LFT(1,planes(i),pin,:))...
    +i*offset),'r','displayname','high flow');
ph(2)=plot(angu,flipud(squeeze(G.LFT(3,planes(i),pin,:)-G.LFT(1,planes(i),pin,:))...
    +i*offset),'g','displayname','mid flow');
ph(3)=plot(angu,flipud(squeeze(G.LFT(4,planes(i),pin,:)-G.LFT(1,planes(i),pin,:))...
    +i*offset),'b','displayname','low flow');
end
legend(ph)
grid on
xlim(G.d_angle/2*[-1,1])
xticks([-90 -45,0,45 90])
xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
%% void faction

fh=figure(8334);clf
cas=[2];
pla=[10]

imshow( -(rec(:,:,cas,pla)-rec(:,:,1,pla))',[])
set(gca,'YDir','normal')
% set(fh,'units','normalized',...
%     'outerposition',[0 0 1 1]);

%fh=figure(3425);clf
P3=f42_ChannelPhanTransform_w_film();

P3=f42_RotPhan(P3,deg2rad(G.cen.prot));
P3=f42_PhanZoom(P3,G.cen.pzoom);
P3=f42_PhanMove(P3,G.cen.px,G.cen.py);
f42_PlotPhan(P3,gca,[1 0 0],[1]);
axis off
% set(gca,'units','pixel',...
%     'Position',[0 0 T.Rec.recsize T.Rec.recsize]);
vm=getframe(gca);

vmask=im2bw(vm.cdata(:,:,2),0.7);
vmask=~imfill(vmask,[1 1]);
vmask=~imfill(vmask,[1 1]);
vmask=~imfill(vmask,[1 1]);
vmask=~imfill(vmask,[1 1]);
vmask=~imfill(vmask,[1 1]);
vmask=~imfill(vmask,[1 1]);
%%
fh2=figure(3426);clf;
imshow(vmask);
title('vmask')
hold on
%f42_PlotPhan(P3,gca,[1 0 0],[1]);
plot(G.cen.cen([1 2 3 4 1],1),G.cen.cen([1 2 3 4 1],2),'k');
plot(G.cen.cen([1 3 2 4],1),G.cen.cen([1 3 2 4],2),'k');
vm2=getframe(gca);

vmask2=im2bw(vm2.cdata);
figure(256);clf;
imshow(vmask2,[]);
title('vmask2')
figure(257);clf;
[L,n] = bwlabel(vmask2);
imshow(L,[])
title('L')

%%
voidmean=zeros(n,G.raw.ncas,G.raw.maxp);


for cas=1:G.raw.ncas
    
    for pla=1:G.raw.maxp
        tom=rec(:,:,cas,pla);
        for reg=1:n
            voidmean(reg,cas,pla)=mean(tom(L==reg));
        end
    end
end

%% plot some void profiles
inner=[3 4 5 6];
outer=[1 2 7 8];
nospac=[1:3,6:46];

cas=3;
reg=6;
figure(4572);clf;
ax=gca()
ylims=[0,1.8]*1e-6;
xlims=[1,44];

subplot(2,2,1)
hold on
for reg=inner
    plot(squeeze(voidmean(reg,2,nospac)-voidmean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end
grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title('inner regions. High flow')
ylabel('entrainment')


subplot(2,2,3)
hold on
for reg=outer
    plot(squeeze(voidmean(reg,2,nospac)-voidmean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end
grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title('outer regions')
ylabel('entrainment')
xlabel('channel length')

subplot(2,2,2)
hold on
for reg=inner
    plot(squeeze(voidmean(reg,3,nospac)-voidmean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end
grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title('inner regions. Low flow')
ylabel('entrainment')


subplot(2,2,4)
hold on
for reg=outer
    plot(squeeze(voidmean(reg,3,nospac)-voidmean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end
grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title('outer regions')
ylabel('entrainment')
xlabel('channel length')

%%
figure(1346);clf
imshow((rec(:,:,cas,pla)-rec(:,:,1,pla))'+1e-5*vmask,[])

void=(rec(:,:,cas,pla)-rec(:,:,1,pla))'.*vmask;
figure(257);clf;
imshow(void,[])
title('void signal')

%% 
figure(4835);clf
ax=cla();
hold on
 cases=[2 3 4]
tstr={'high flow','mid flow','low flow'} 
for i=1:3 
subplot(1,3,i)
hold on
for reg=inner
    plot(squeeze(voidmean(reg,cases(i),nospac)-voidmean(reg,1,nospac)),...
        'displayname',sprintf('reg. %d',reg));
end

grid on
ylim(ylims)
xlim(xlims)
legend('location','best')
title(tstr{i})
ylabel('entrainment')
end
%% de-centering plot


cas=8;
rep=1;
fname=sprintf('%s%02d_%02d.mat','1_corr_',cas,rep);
fpath=sprintf('%s%s',T.d.DataPath,fname);
d=G.BoLoad(fpath,T);

% crop z
d=d(:,:,T.q360.startframe:(...
    T.q360.startframe+T.q360.nFrames(1,1)-1));

% image correction: pixels & filter
corr=G.ImCorrFilt(mean(d,3),T);
% apply image correction
d=G.ApplyImCorr(d,corr,T);

% crop x & y
d=d(T.Cen.range,:,:);

% quadrant-rotate, fliplr, log
d=flip(-log(circshift(d,T.Cen.quadrot,3)),3);


%%
pla=600;
[ax,fh]=pub.decenter(T,d)

%% center shift  scan

sino=squeeze(d(:,pla,:));

% center exact
%baseshift(plane)
nshift=41; % number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);

mask=G.MakeTomoMask('fan',G.fraccircshift(...
                sino,extrashift(21)),recsize,ang,detPitch,src,det,...
                MaskThresh,0.000005,[150 150]);
            
            
% reconstruct
rec=zeros(recsize,recsize,nshift);
for i =1:nshift
rec(:,:,i)=a.FBPexplFan(...
    G.fraccircshift(sino,extrashift(i))',...
    recsize,ang,detPitch,src,det);
    [q(i),~]=qc.VarianceQuality(rec(:,:,i),...
        mask);
end
%%
figure(23451);clf
imshow(mask,[])
title('mask')
fh=figure(2314);clf
mask2=mask;
mask2(round(recsize/2):end,1:round(recsize/2))=0;
imshowpair(imgradient(rec(:,:,21))',mask2'>1);
set(gca,'YDir','normal')
text(153,250,'image gradient','Color',[0 1 0],'FontSize',12)
text(220,230,'mask','Color',[1 0.2 1],'FontSize',12)
ax=gca;
ax.Units='pixel';
ax.Position=[0 0 recsize recsize];
fh.Units='pixel';
fh.Position=[50 50 recsize recsize];
fh.Units='centimeters';
fh.PaperPosition=fh.Position;
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperUnits', 'centimeters')
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

fname=sprintf('Qualitymask1');
print(sprintf('%s%s.pdf',T.G.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.G.saveto,fname))
%%
fh=figure(4894);clf;
plot(extrashift,q,'x')
grid on
ylh=xlabel('shift [pixel]')
xlh=ylabel('$\sigma(\nabla Im\mid_{Mask})$','Interpreter','latex')
pub.BFxlab(xlh,T.F)
pub.BFylab(ylh,T.F)


fh.Units='centimeters';
fh.Position=[2 2 7.1173 7.1173];

set(fh, 'PaperUnits', 'centimeters')
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

ax = gca;
pub.BFaxis(ax,T.F)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fname=sprintf('Qualitymask2');
print(sprintf('%s%s.pdf',T.G.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.G.saveto,fname))
%%
fh=figure(48941);clf;
plot( T.Cen.eshift(:,cas,rep))
hold on
plot( pla,T.Cen.eshift(pla,cas,rep),'x')

grid on
ylh=xlabel('z-plane')
xlh=ylabel('best shift [pixel]')
xlim([0 1024])
ylim(1.2*[-1,1])
pub.BFxlab(xlh,T.F)
pub.BFylab(ylh,T.F)


fh.Units='centimeters';
fh.Position=[2 2 7.1173 7.1173];

set(fh, 'PaperUnits', 'centimeters')
fh.PaperSize=fh.Position(3:4);
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

ax = gca;
pub.BFaxis(ax,T.F)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];
fname=sprintf('Qualitymask3');
print(sprintf('%s%s.pdf',T.G.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.G.saveto,fname))

%%
fh=figure(49846);clf;
im=[];
xrange=30:80;
yrange=50:120;
c1=[];
for i=1:5:nshift
   im=cat(1,im,rec(xrange,yrange,i)); 
   c1=cat(1,c1,i);
end

im2=[];
c2=[];
for i=17:25
   im2=cat(1,im2,rec(xrange,yrange,i)); 
   c2=cat(1,c2,i);
end

im=cat(2,im2,im);
imshow(im',[])
set(gca,'YDir','normal')

for i=1:9
    text(25+(i-1)*51,95,sprintf('%2.1f',extrashift(c1(i))),...
        'color',[0.99 1 1],'FontSize',12,'FontName', 'Helvetica')
    text(25+(i-1)*51,25,sprintf('%2.1f',extrashift(c2(i))),...
        'color',[.99 1 1],'FontSize',12,'FontName', 'Helvetica')
end

ax=gca;
ax.Units='pixel';
ax.Position(1:2)=[0 0];

fh.Units='pixel';
fh.Position=[50 50 ax.Position(3:4)];
fh.Units='centimeters';


set(fh, 'PaperUnits', 'centimeters')
set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
fh.PaperSize=fh.Position(3:4);

    fname=sprintf('Qualitymask4');
print(sprintf('%s%s.pdf',T.G.saveto,fname),'-dpdf')
open(sprintf('%s%s.pdf',T.G.saveto,fname))