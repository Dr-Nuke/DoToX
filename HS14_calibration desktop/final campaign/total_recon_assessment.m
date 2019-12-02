% this file will doe the final assessment of the algos and provide
% publication grade figures.

% topics to cover:
% FBP vs SART vs SIRT
% parallel vs fan vs cone
% 
% make mask
% explain quality criterion
% -centering search
% -runtime vs iteration vs result
% 

%%
if exist('T','var')
    clearvars -except T
    
else
    clear
    load('T')
end

clc
close all
format compact

% add astra paths
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')

extracrop=20;
cas=[1,6,8]; % what case {1-10}
cas2=[1,2,4];% what case (recon format {1-4})
rep=[1,9]; % what repetitions
pla=[50 100 ,250, 310, 350, 512, 640, 1000]; %sino planes
planam={'pins','before spacer','spacer grid','vanes','above spacer',...
    'mid plane','640', 'top'};
ang=deg2rad(T.Rec.angles); % projection angles list
xrange=T.Cen.start+extracrop:T.Cen.stop-extracrop; % rough cropping around channel
xrange=1:269;
quadrot=70; % rotate to have the channel quadrants aligned
sinorange=T.q360.startframe:T.q360.startframe+T.q360.NFGuess-1; % 360°crop
sinorange=1:1357;
recsize=length(xrange);


tra_raw=[];
for c=1:length(cas)
    for r=1:length(rep)
        d=f.BoLoad(T.fnames.corr{cas(c),rep(r)},T);
        tra_raw(:,:,:,c,r)=d(xrange,pla,sinorange);
        shift_raw(:,c,r)=T.Cen.fitshift(cas2(c),pla);
        
    end
end
save('tra_raw','tra_raw');
save('shift_raw','shift_raw')
clear d

% cheks:
load('tra_raw')
load('shift_raw')
size(tra_raw)
nframes=size(tra_raw,3);
tra=flip(-log(circshift(tra_raw,quadrot,3)),3);

%make some easy sinos
%%
r=1;

detPitch=0.127; % detector pixel spacing
src=950; % source-axis distance
det=50; % detector-axis distance

for c=1:length(cas)
    for p=1:length(pla)
        check1_sino(:,:,c,p)=f.fraccircshift(squeeze(tra(:,p,:,c,1)),...
            -shift_raw(p,c,r));
        check1_rec(:,:,c,p,1)=a.FBPexpl(check1_sino(:,:,c,p)',recsize,ang);
        check1_rec(:,:,c,p,2)=a.FBPexplFan(check1_sino(:,:,c,p)',...
                    recsize,ang,detPitch,src,det);
        
    end
end
%% plot

c=1

for p=1:length(pla)
figure(100+p);clf
im1=check1_rec(:,:,c,p,1);
im1max=mean(im1(:));
%im2=ones(recsize,recsize);
im2=check1_rec(:,:,c,p,2);
im2max=mean(im2(:));

imshow(cat(1,im2/im2max,im1/im1max),[])
title(sprintf('c%d r%d p%d top par bot fan ',...
    c,r,pla(p)))
set(gca,'YDir','normal')
    
end


%%

% make the mask for each algo
nmeth=3; % parallel fan and cone
mask=zeros(recsize,recsize,nmeth);
recmask=zeros(recsize,recsize,nmeth);

% define the reference case:
c=1; % case
p=2; % plane
r=1; % rep
thresh=0.003;

% center roughly
mask_sino=f.fraccircshift(squeeze(tra(:,p,:,c,1)),...
            -shift_raw(p,c,r));
mask=[];
recmask(:,:,1)=a.FBPexpl(mask_sino',recsize,ang);
[~,mask(:,:,1),gm1(:,:,1)]=qc.kern(recmask(:,:,1),thresh,[150,150]);

recmask(:,:,2)=a.FBPexplFan(mask_sino',recsize,ang,detPitch,src,det);
[~,mask(:,:,2),gm1(:,:,2)]=qc.kern(recmask(:,:,2),0.0000038,[150,150]);
figure(34);clf; imshowpair(mask(:,:,2),mask(:,:,1))
title(sprintf('two masks for para & fan %d %d',...
    sum(sum(mask(:,:,1))),sum(sum(mask(:,:,2)))))
%%
pub.QualityMask(recmask,gm1,mask,T)
%%


%% n- step centering
close all
p=7;
c=2;
r=1;
nmeth=2; % parallel & fan beam
nshift=21; % number of tests
minshift=-0.5; % in pixels
maxshift=0; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);

nsc_sinogram=[];
nsc_rec=[];
nsc_v=[];


for i=1:nshift
%     nsc_sinogram(:,:,i)=f.fraccircshift(cat(1,zeros(10,1357),sg,zeros(10,1357)),...
%         -extrashift(i));
    nsc_sinogram(:,:,i)=f.fraccircshift(squeeze(tra(:,p,:,c,r)),...
            -(shift_raw(p,c,r)+extrashift(i)));
    m=1; % parallel    
    nsc_rec(:,:,i,m)=a.FBPexpl(nsc_sinogram(:,:,i)',recsize,ang);
    [nsc_v(i,m),~]=qc.VarianceQuality(nsc_rec(:,:,i,m),mask(:,:,m));
    
    m=2; % fan beam
    nsc_rec(:,:,i,m)=a.FBPexplFan(nsc_sinogram(:,:,i)',recsize,ang,...
        detPitch,src,det);
    [nsc_v(i,m),~]=qc.VarianceQuality(nsc_rec(:,:,i,m),mask(:,:,m));
end

% plots
%close all
%xr=1:recsize;
xr=20:120;
%yr=1:recsize;
yr=20:120;
nsctotim=[];
im4=[];
xstr= [];
figure(301);clf
for i=1:nshift
    im1=nsc_rec(xr,yr,i,1);
    maxim1=max(im1(:));
    
    im2=nsc_rec(xr,yr,i,2);
    maxim2=max(im2(:));
    
    im3=(cat(1,im2'/maxim2,im1'/maxim1));
    im4=cat(2,im4,im3);
    xstr=sprintf('%s %.2f',xstr,extrashift(i));
end
imshow(im4,[]);
set(gca,'YDir','normal')
title(sprintf('c%d r%d p%d top parallel bottom fan beam, shift from %.2f to %.2f',...
    c,r,pla(p),minshift,maxshift))
xlabel(xstr)


[vmax,vind]=max(nsc_v(:,1));
figure(2356);clf
imshowpair(mask(:,:,1),imgradient(nsc_rec(:,:,vind,1)))
title(sprintf('%d imshowpair parallel',vind))

[vmax,vind]=max(nsc_v(:,2));
figure(2357);clf
imshowpair(mask(:,:,2),imgradient(nsc_rec(:,:,vind,2)))
title(sprintf('%d imshowpair fan',vind))

figure(34);clf
yyaxis left
plot(extrashift,nsc_v(:,1),'xr','Displayname','parallel')

yyaxis right
plot(extrashift,nsc_v(:,2),'+b','Displayname','fan')
grid on
legend()




%% make one quality plot for parallel vs fan beam


%% test fibonacci centering vs classic

% ace parameter etc
c=2; % case
r=1; % rep
p=5; % plane
m=1; % method i: even distributed points (vs. fibonacci)
al=1; % alog 1: fbp

% more parameters
nshift=50; % number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift); % array of extra shifts
nalgo=3; % FBP SIRT SART 
nmeth=2; % even or fibonacci


%fibonacci stuff
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
epsilon=0.01;               % accuracy value // intervall size lower bound
                  

% preallocation
v=zeros(nshift,length(cas),length(rep),length(pla),nmeth,nalgo);% variance quality indicator
xlist=zeros(nshift,length(cas),length(rep),length(pla),nalgo);% variance quality indicator
t=zeros(length(cas),length(rep),length(pla),nmeth,nalgo); % time requirement

rec=zeros(recsize,recsize,nshift,length(cas),length(rep),...
    length(pla),nmeth,nalgo,'single');% recon container
x1n=@(a,b,tau)a+(1-tau)*(b-a);   % lefz mid point
x2n=@(a,b,tau)a+tau*(b-a);      % riht mid point

fprintf('check regular vs fibonacci spacing....')
for c=1:2 % cases
    for r=1:2%2 % reps
        for p=[2 6 8] % planes
            disp([c,r,p]) % debug
            sino=squeeze(tra(:,p,:,c,r));
            shift=(T.Cen.fitshift(cas2(c),pla(p))); % original calculatd shift
            
            % regular spacing
            m=1; % mode 1: recular
            tic
            for i=1:nshift
                sinogram=f.fraccircshift(sino,...
                    -(shift+extrashift(i)));
                rec(:,:,i,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                    recsize,ang,detPitch,src,det);
                [v(i,c,r,p,m,al),~]=qc.VarianceQuality(rec(:,:,i,c,r,p,m,al),mask(:,:,2));
            end

            t(c,r,p,m,al)=toc;
            
            % fibonacci spacing
            tic
            m=2; % mode 2: fibonacci
            aa=minshift; % intervall bounds
            b=maxshift;
            i=2;  % number of (initial) iterations, counter
            
            x1=x1n(aa,b,tau);             % computing initial mid x values
            x2=x2n(aa,b,tau);
            xlist(1:2,c,r,p,al)=[x1,x2]; % list for plotting
            
            % find the two initial values
            %1
            sinogram=f.fraccircshift(sino,...
                -(shift+x1));
            rec(:,:,1,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                recsize,ang,detPitch,src,det);;
            [f1,~]=qc.VarianceQuality(rec(:,:,1,c,r,p,m,al),mask(:,:,2));
            
            %2
            sinogram=f.fraccircshift(sino,...
                -(shift+x2));
            rec(:,:,2,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                recsize,ang,detPitch,src,det);
            [f2,~]=qc.VarianceQuality(rec(:,:,2,c,r,p,m,al),mask(:,:,2));
   
            v(1:2,c,r,p,m,al)=[f1,f2];
            while ((abs(b-aa)>epsilon) && (i<nshift))
                i=i+1;
                if(f1>f2)
                    b=x2;
                    x2=x1;
                    x1=x1n(aa,b,tau);
                    xlist(i,c,r,p,al)=x1;
                    f2=f1;
                    
                    sinogram=f.fraccircshift(sino,...
                        -(shift+x1));
                    rec(:,:,i,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                        recsize,ang,detPitch,src,det);
                    [f1,~]=qc.VarianceQuality(rec(:,:,i,c,r,p,m,al),mask(:,:,2));
                    v(i,c,r,p,m,al)=f1;
                else
                    aa=x1;
                    x1=x2;
                    x2=x2n(aa,b,tau);
                    xlist(i,c,r,p,al)=x2;
                    
                    f1=f2;
                    sinogram=f.fraccircshift(sino,...
                        -(shift+x2));
                    rec(:,:,i,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                        recsize,ang,detPitch,src,det);
                    [f2,~]=qc.VarianceQuality(rec(:,:,i,c,r,p,m,al),mask(:,:,2));
                    v(i,c,r,p,m,al)=f2;

                end
                
                
            end
            t(c,r,p,m,al)=toc;

        end
    end
end

%%

for c=1:2 % cases
    for r=1:2 % reps
        for p=[2 6 8] % planes
            figure();clf;
            plot(extrashift,v(:,c,r,p,1,al),'xr','displayname',sprintf('c%d r%d p%d reg',c,r,p))
            hold on
            grid on
            plot(xlist(:,c,r,p,al),v(:,c,r,p,2,al),'xb','displayname',sprintf('c%d r%d p%d fib',c,r,p))
            legend()
            title(sprintf('c%d r%d p%d',c,r,p))
        end
    end
end
%%
repcat=[]
for i =1:5:nshift
    repcat=cat(2,repcat,rec(50:260,50:260,i,1,1,2,1,1));
    
end
figure(451);clf;
imshow(repcat,[])

% hold on
% plot(extrashift,v(:,2),'+b','displayname',sprintf('SIRT, %.1fs, 100 iter.',t(2) ))
% plot(extrashift,v(:,3),'*g','displayname',sprintf('SART, %.1fs, 1800 iter.',t(3) ))

% xlabel('centering shift')
% ylabel('gradient sum of squares')
% grid on
% title('optimization for centering')
% legend()


%% paralles vs fan beam
% set parameters
c=1; %case
r=1; % rep
p=2;    % plane
al=1;   % algo


% get best centering shift
[vmax,vind]=max(v(:,c,r,p,2,al));
shift=(T.Cen.fitshift(cas2(c),pla(p)))+xlist(vind,c,r,p,al); 

                
%preallocate                
pvf_rec=zeros(recsize,recsize,length(pla),2);
pvf_v=zeros(length(pla),2);


for p=1:length(pla)   
    sinogram=f.fraccircshift(squeeze(tra_cor(:,p,:,c,r)),-shift);   
    pvf_rec(:,:,p,1)=a.FBPexpl(sinogram',recsize,ang); 
    pvf_rec(:,:,p,2)=a.FBPexplFan(sinogram',recsize,ang,detPitch,src,det);
    [pvf_v(p,1),~]=qc.VarianceQuality(pvf_rec(:,:,p,1),mask);
    [pvf_v(p,2),~]=qc.VarianceQuality(pvf_rec(:,:,p,2),mask);

end
         
% plot
%
xrange=50:130;
yrange=190:270;
pvf_plotrec=[];
for p=1:length(pla)    
   
    pvf_rectm(p,1)=mean(mean(pvf_rec(:,:,p,1)));
    pvf_rectm(p,2)=mean(mean(pvf_rec(:,:,p,2)));
    pvf_plotrect=cat(1,pvf_rec(xrange,yrange,p,1)'/pvf_rectm(p,1),...
                     pvf_rec(xrange,yrange,p,2)'/pvf_rectm(p,2));
    pvf_plotrec=cat(2,pvf_plotrec,pvf_plotrect);
    
end
%
figure(123);clf;
imshow(pvf_plotrec,[]); set(gca,'YDir','normal');
figure(124);clf;


%legend('parallel','fan beam')
xlabel('planes 50 100 350 640 1000')
ylabel('quality criterion')


%% new parallel vs. fanbeam vs cone beam via equidistanr
c=2;
r=2;

c=2;
r=2;
p=4;
nshift=31; % max number of tests
minshift=-2; % in pixels
maxshift=1; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift); % array of extra shifts
nmeth=3; % parallel, fan, cone
nalgo=3;

% preallocation
pfc_v=zeros(nshift,length(pla),nalgo);% variance quality indicator
pfc_rec=zeros(recsize,recsize,nshift,length(pla),nalgo,'single'); %recons
pfc_sinogram=zeros(size(sino,1),size(sino,2),nshift,length(pla),'single');


sino=squeeze(tra_cor(:,p,:,c,r));
shift=(T.Cen.fitshift(cas2(c),pla(p))); % original calculatd shift

for i=1:nshift
    disp(i)
    pfc_sinogram(:,:,i,p)=f.fraccircshift(sino,...
        -(shift+extrashift(i)));
    
    % 1) fbp
    al=1;
    pfc_rec(:,:,i,p,al)=a.FBPexpl(pfc_sinogram(:,:,i,p)',...
        recsize,ang);
    [pfc_v(i,p,al),~]=qc.VarianceQuality(pfc_rec(:,:,i,p,al),mask);
    
    % 1) fan
    al=2;
    pfc_rec(:,:,i,p,al)=a.FBPexplFan(pfc_sinogram(:,:,i,p)',...
        recsize,ang,detPitch,src,det);
    [pfc_v(i,p,al),~]=qc.VarianceQuality(pfc_rec(:,:,i,p,al),mask);

 
    
end
%% plot
pfc_imrec=[];
tstr='';
% for i=1:5:nshift
%     
%     pfc_rectm(i,1)=max(max(pfc_rec(:,:,i,p,1)));
%     pfc_rectm(i,2)=max(max(pfc_rec(:,:,i,p,2)));
%     
%     pfc_imrect=cat(1,pfc_rec(:,:,i,p,2)'/pfc_rectm(i,2),...
%                      pfc_rec(:,:,i,p,1)'/pfc_rectm(i,1));
%     pfc_imrec=cat(2,pfc_imrec,pfc_imrect);
%     tstr=sprintf('%s%2d ',tstr,i);
% end
close all
for i=1:nshift
    fh=figure();
    fh.Position=[1007 535 560 420];
    imshow(pfc_rec(:,:,i,p,1),[])
    title(sprintf('par %2d, %0.3f',i,pfc_v(i,p,1)))
    
    fh=figure();
    fh.Position=[422 547 560 420];
    imshow(pfc_rec(:,:,i,p,1),[])
    title(sprintf('fan %2d, %0.3f',i,pfc_v(i,p,1)))
end
    




%% new parallel vs. fanbeam vs cone beam via golden cut method
% pfc
c=2;
r=2;
nshift=30; % max number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift); % array of extra shifts
nmeth=3; % parallel, fan, cone

%fibonacci stuff
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
epsilon=0.01;               % accuracy value // intervall size lower bound

% preallocation
pfc_v=zeros(nshift,length(pla),nmeth);% variance quality indicator
pfc_xlist=zeros(nshift,length(pla),nmeth);% variance quality indicator
pfc_rec=zeros(recsize,recsize,nshift,length(pla),nmeth,'single'); %recons
pfc_t=zeros(length(pla),nmeth); % time requirement

for p=1:length(pla)
    % chose plane sino
    sino=squeeze(tra_cor(:,p,:,c,r));
    shift=(T.Cen.fitshift(cas2(c),pla(p))); % original calculatd shift
            
    
    % 1 parallel
    m=1
    aa=minshift; % intervall bounds
    b=maxshift;
    i=2;  % number of (initial) iterations, counter
    
    x1=x1n(aa,b,tau);             % computing initial mid x values
    x2=x2n(aa,b,tau);
    pfc_xlist(1:2,p,m)=[x1,x2]; % list for plotting
   
    sinogram=f.fraccircshift(sino,...
        -(shift+x1));
    pfc_rec(:,:,1,p,m)=a.FBPexplFan(sinogram',...
        recsize,ang,detPitch,src,det);
    [f2,~]=qc.VarianceQuality(pfc_rec(:,:,1,p,m),mask);
    
    v(1:2,c,r,p,m,al)=[f1,f2];
    while ((abs(b-aa)>epsilon) && (i<nshift))
        
        i=i+1;
        if(f1>f2)
            b=x2;
            x2=x1;
            x1=x1n(aa,b,tau);
            xlist(i,c,r,p,al)=x1;
            f2=f1;
            
            sinogram=f.fraccircshift(sino,...
                -(shift+x1));
            rec(:,:,i,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                recsize,ang,detPitch,src,det);
            [f1,~]=qc.VarianceQuality(rec(:,:,i,c,r,p,m,al),mask);
            v(i,c,r,p,m,al)=f1;
        else
            aa=x1;
            x1=x2;
            x2=x2n(aa,b,tau);
            xlist(i,c,r,p,al)=x2;
            
            f1=f2;
            sinogram=f.fraccircshift(sino,...
                -(shift+x2));
            rec(:,:,i,c,r,p,m,al)=a.FBPexplFan(sinogram',...
                recsize,ang,detPitch,src,det);
            [f2,~]=qc.VarianceQuality(rec(:,:,i,c,r,p,m,al),mask);
            v(i,c,r,p,m,al)=f2;
            
        end
        
        
    end
end

     
%% speed test
p=4;
c=1;
r=2;
sino=squeeze(tra_cor(:,p,:,c,r));

fprintf('parallel...')
tic

for i = 1:nshift
    st_rec=a.FBPexpl(f.fraccircshift(sino,...
        -(shift+extrashift(i)))', recsize,ang); 
    
end
toc
  
fprintf('fan beam...')
tic
for i = 1:nshift

    st_rec=a.FBPexplFan(f.fraccircshift(sino,...
        -(shift+extrashift(i)))',recsize,ang,detPitch,src,det);
    
end
toc


sino=squeeze(tra_cor(:,p,:,c,r))

fprintf('GPU parallel...')
tic

for i = 1:nshift
    st_rec=a.FBPexpl(f.fraccircshift(sino,...
        -(shift+extrashift(i)))', recsize,ang); 
    
end
toc
  
fprintf('GPU fan beam...')
tic
for i = 1:nshift

    st_rec=a.FBPexplFan(f.fraccircshift(sino,...
        -(shift+extrashift(i)))',recsize,ang,detPitch,src,det);
    
end
toc





%% histogram comparison







            
