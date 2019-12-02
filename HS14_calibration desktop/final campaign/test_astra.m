% test the astra functionality
% do the FBP algo and then adjust the iteratives until they match the time
% consumtion



% clear 
% close all
% % first load some data
% load('T')

addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\mex')
addpath('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\tools')
d=f.BoLoad(T.fnames.corr{6,1},T);
d=d(:,:,T.q360.startframe:T.q360.startframe+T.q360.NFGuess-1);
%

%create mask
maskplane = 100;
height=1;
thresh=0.002;
ang=T.Rec.angles;
% log, squeez, and mean of stack height
sinomask=-log(squeeze(mean(d(173:491,maskplane:maskplane+height-1,:),2)));


% shift to have the quadrants aligned
sinomask=circshift(sinomask,T.Rec.SinoRotShift,2);
%create the unique mask
sinomask=f.fraccircshift(sinomask,-(T.Cen.fitshift(2,maskplane)));
recmask=iradon(sinomask,ang,'spline','Hann',1,319);
[~,mask,gm1]=qc.kern(recmask,thresh,[160,160]);

figure(55);clf;
imshow(cat(1,mask',gm1'/max(gm1(:)),recmask'/max(recmask(:))),[])
set(gca,'YDir','normal')
title('recon - grad - mask')

A=a.CreateAstraDataStruct();
ang=(A.ang);
%% centering optimization


plane=521;
nshift=12; % number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);
A=a.CreateAstraDataStruct();
ang=(A.ang);

tic
sino=(-log(squeeze(d(173:491,640,:))));
sino=circshift(sino,T.Rec.SinoRotShift,2); % the quadrant-align-rotation
v=zeros(size(3,nshift));

shift=(T.Cen.fitshift(2,640));
fprintf('FBP: ')
q=1;
tic
for i=1:nshift
    sinogram=f.fraccircshift(sino,-(shift+extrashift(i)))';
    rec(:,:,i,q)=a.FBPexpl(sinogram,319,ang);
    [v(i,q),~]=qc.VarianceQuality(rec(:,:,i,q),mask);
end
t(1)=toc;
fprintf('%.2f\n',t(1));

fprintf('SIRT: ')
q=2;
tic
for i=1:nshift 
    sinogram=f.fraccircshift(sino,-(shift+extrashift(i)))';
    rec(:,:,i,q)=a.SIRTexpl(sinogram,319,ang,100);
    [v(i,q),~]=qc.VarianceQuality(rec(:,:,i,q),mask);
    
end
t(2)=toc;
fprintf('%.2f \n',t(2));

fprintf('SART: ')
q=3;
tic
for i=1:nshift 
    sinogram=f.fraccircshift(sino,-(shift+extrashift(i)))';
    rec(:,:,i,q)=a.SARTexpl(sinogram,319,ang,1800);
    [v(i,q),~]=qc.VarianceQuality(rec(:,:,i,q),mask);
end
t(3)=toc;
fprintf('%.2f\n',t(3));

disprec=[];
for i = 1:2
    dr=[];
    for j=1:6
        k=j+6*(i-1);
        dr=cat(2,dr,rec(:,:,k));
        
    end
    disprec=cat(1,disprec,dr);
end
figure(59);clf
imshow(imgradient(disprec),[])




figure(60);clf
plot(extrashift,v(:,1),'xr','displayname',sprintf('FBP, %.1fs',t(1) ))
hold on
plot(extrashift,v(:,2),'+b','displayname',sprintf('SIRT, %.1fs, 100 iter.',t(2) ))
plot(extrashift,v(:,3),'*g','displayname',sprintf('SART, %.1fs, 1800 iter.',t(3) ))

xlabel('centering shift')
ylabel('gradient sum of squares')
grid on
title('optimization for centering')
legend()


[vmax,vind]=max(v);

figure(61);clf

drec=[];
for i = 1:3
    rect=rec(:,:,vind(i),i);
    recd=imgradient(rect);
    recmax=max(rect(:));
    recdmax=max(recd(:));
    
    drect=cat(2,rect/recmax,recd/recdmax);
    drec=cat(1,drec,drect);
end

imshow(drec,[])
title('top:FBP middle: SIRT botom: SART')
% create properties
%A=a.CreateAstraDataStruct(T);

%% check fr iteration dependence

iters=[10 20 50  100 150 200 250 300 350 400 450 500 1000 2000];
niter=length(iters);
sino1=f.fraccircshift(sino,-(shift+extrashift(vind(2))))'; % SIRT optimal
sino2=f.fraccircshift(sino,-(shift+extrashift(vind(3))))'; % SART optimal
iv=[];
it=[];
ylabstr=[]
irec=[]
for i=1:niter
    disp(i);
    
    tic
    irec(:,:,i,2)=a.SARTexpl(sino2,319,ang,iters(i));
    it(i,2)=toc;
    
    tic
    irec(:,:,i,1)=a.SIRTexpl(sino1,319,ang,iters(i));
    it(i,1)=toc;
    

    
    [iv(i,1),~]=qc.VarianceQuality(irec(:,:,i,1),mask);
    [iv(i,2),~]=qc.VarianceQuality(irec(:,:,i,2),mask);
    ylabstr=[ylabstr sprintf('%d ',iters(i))];
    
end
    
    
figure(63);clf;
yyaxis left
plot(iters,iv(:,1),'+b','displayname',sprintf('SIRT quality' ))
hold on
plot(iters,iv(:,2),'*g','displayname',sprintf('SART quality' ))
ylabel('gradient sum of squares')

yyaxis right
plot(iters,it(:,1),'ob','displayname',sprintf('SIRT gpu time' ))
plot(iters,it(:,2),'sg','displayname',sprintf('SIRT gpu time' ))

set(gca, 'XScale', 'log')
xlabel('# iterations')
ylabel('computation time [s]')
grid on
title('optimization for iterations & Algo')
legend('location','northwest')

figure(64);clf

drec=[];
for i = 1:niter
    drect=[];
    rect=irec(:,:,i,1);
    recd=imgradient(rect);
    recmax=max(rect(:));
    recdmax=max(recd(:));
    drect=cat(2,rect/recmax,recd/recdmax);
    
    rect=irec(:,:,i,2);
    recd=imgradient(rect);
    recmax=max(rect(:));
    recdmax=max(recd(:));
    
    drect=cat(2,drect,rect/recmax,recd/recdmax);
    
    
    drec=cat(1,drec,drect);
end

imshow(drec,[])
title(sprintf('from left to right: \n SIRT, grad(SIRT), SART, grad(SART)'))
ylabel(['iterations: ' ylabstr])
set(get(gca,'YLabel'),'Rotation',270)

%% SART vs FBP
imrecdif=[];
dif=[]
recdif=[];%zeros(size(irec));
for i=1:niter
    recdif(:,:,i,1)=rec(:,:,vind(1),1)-irec(:,:,i,1); %sirt
    recdif(:,:,i,2)=rec(:,:,vind(1),1)-irec(:,:,i,2); %sart
    
    rect=cat(2,rec(:,:,vind(1),1)',irec(:,:,i,1)',recdif(:,:,i,1)',irec(:,:,i,2)',recdif(:,:,i,2)');
    imrecdif=cat(1,imrecdif,rect);
    dif(i,1)=sum(sum(recdif(:,:,i,1).^2)); %Sirt
    dif(i,2)=sum(sum(recdif(:,:,i,2).^2)); %Sart
end
figure(67);clf;
imshow(imrecdif,[]);
xlabel('FBP (ref), SIRT, SIRT diff, SART, SART diff')
ylabel(['iterations: ' ylabstr])
set(get(gca,'YLabel'),'Rotation',270)

figure(68);clf;
plot(iters,dif(:,1),'+b','displayname',sprintf('SIRT' ))
hold on
plot(iters,dif(:,2),'+g','displayname',sprintf('SART' ))
grid on
legend()
ylabel('squared deviation from FBP')
xlabel('iterations')
title('SIRT and SART iteration study')

set(gca, 'XScale', 'log')
%set(gca,'YDir','normal');


%% parallel vs fan beam
%sino1=f.fraccircshift(sino,-(shift+extrashift(vind(2))));
res=0.127;
src=950;
det=50;
fanrec=a.FBPexplFan(sino1,319,ang,0.127,950,50);
figure(75);clf
%imshow(cat(2,fanrec,rec(:,:,vind(1),1)),[])
imshow(fanrec,[])
title('fan beam')
colorbar

figure(76);clf
imshow(rec(:,:,vind(1),1),[])
colorbar
title('parallel')


%% create 4 slices for micha

plane=512;
planes=512:521;
cases= [1 6 8]; % case for preprocessing
fnamelabel={'NoFilm' 'BigFilm' 'SmallFilm'}
cases2=[1 2 4]; % for reconstructing
reps=[1 2 3 4 5 6 7 8 9]; % repetition
Michasino=[];
xrange=193:461;
Mrecsize=length(xrange);
if 0
for cc=1:length(cases)
    for r=reps
   % read in file
    d=f.BoLoad(T.fnames.corr{cases(cc),reps(r)},T);
     
   %extract sino
   Michasino(:,:,:,cases2(cc),r)=d(xrange,plane:plane+9,...
       T.q360.startframe:T.q360.startframe+T.q360.NFGuess-1);
    end
end
save('Michasino','Michasino')
else
    load('Michasino')
end

rawsino2=flip(-log(circshift(Michasino,quadrot,3)),3);
%%


MMask=f.MakeTomoMaskPar(f.fraccircshift(squeeze(rawsino2(:,1,:,1,1)),...
    -T.Cen.fitshift(cases2(cc),planes(p))),Mrecsize,ang,0.003,[150 150]);
semimask=MMask;
semimask(1:round(recsize/2),:)=0;
figure(123);clf;
imshow(MMask,[]);
title('Mask for micha')

imax=20;
xmin=-2;
xmax=2;

Mshift=zeros(length(planes),length(cases),length(reps));
xlist=zeros(imax,length(planes),length(cases),length(reps));
v=zeros(imax,length(planes),length(cases),length(reps));
fnamelabel={'NoFilm' 'BigFilm' 'SmallFilm'};
for cc=1:1%length(cases)
    for r=1:1%length(reps)
        for pp=1:1%length(planes);
            sino=squeeze(rawsino2(:,pp,:,cases2(cc),r));
            baseshift=T.Cen.fitshift(cases2(cc),planes(p));
            [Mshift(pp,cc,r),xlist(:,pp,cc,r),v(:,pp,cc,r)]=f.fibo(...
                sino,baseshift,...
                xmin,xmax,imax,0.007,1,ang,Mrecsize,detPitch,src,det,MMask);
            
            Mrec(:,:,pp,cc,r)=a.FBPexplFan(f.fraccircshift(...
                sino,-Mshift(pp,cc,r))',Mrecsize,ang,detPitch,src,det);
            
            fig=figure(123);clf;
            fig.Position=[100 162 400 800];
            subplot('Position',[0.05,0.5,0.9,0.45]);
            imshow(Mrec(:,:,pp,cc,r),[]);
            title(sprintf('c%d r%d p%03d',cases(cc),r,planes(pp)));
            subplot('Position',[0.1,0.1,0.8,0.3]);
           
            pv=v(:,pp,cc,r);
            pind=~isnan(pv);
            plot(xlist(pind,pp,cc,r),v(pind,pp,cc,r),'xr')
            hold on
            [vmax,vind]=max(v(:,pp,cc,r));
            
            plot(xlist(vind,pp,cc,r),v(vind,pp,cc,r),'+b','Displayname',...
                sprintf('%.3f',xlist(vind,pp,cc,r)))
            title('shift finding')
            grid on
            xlim([xmin,xmax]);
            legend()
            
            fpath='T:\cbolesch\Sino_fuer_micha\';
            fname=sprintf('%sMicha_c%d_r%d_p%03d_%s',fpath,cases(cc),r,planes(pp),fnamelabel{cc});
            print(sprintf('%s',fname),'-dpng');
            savefig(fig,fname)
            
            csvwrite(sprintf('%s.csv',fname),Mrec(:,:,pp,cc,r));
            dlmwrite(sprintf('%s_SpaceSeparated.csv',fname),Mrec(:,:,pp,cc,r),' ')
            
            
            
            
        end
    end
end


    %%

nshift=21; % number of tests
minshift=-2; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);



for i =1:length(cases)   
   % scan for correct centering
   shift=T.Cen.fitshift(cases2(i),plane);
   for j=1:nshift
    sinogram=f.fraccircshift(rawsino2(:,:,i),-(shift+extrashift(j)))';
    recm(:,:,i,j)=a.FBPexpl(sinogram,319,ang);
    [vm(i,j),~]=qc.VarianceQuality(recm(:,:,i,j),mask);
   end
end

%
figure(80);clf;
ax=gca;
hold on;
col=hsv(length(cases));
markes={'x' 'x' '+' '+' 'o' 'o'};

for i =1:length(cases)
    plot(extrashift,vm(i,:),'color',col(i,:),'marker',markes{i},...
        'displayname',sprintf('cas %d rep %d',cases(i),reps(i)))
end   

xlabel('centering shift')
ylabel('gradient sum of squares')
grid on
title('optimization for centering')
legend()
%

[vmaxm,vindm]=max(vm,[],2)

fnamelabel={'NoFilm' 'NoFilm' 'BigFilm' 'BigFilm' 'SmallFilm' 'SmallFilm'}

for i =1:length(cases)
    fname=sprintf('recon_cas%d_rep%d_plane%d_%s.csv',cases(i),reps(i),plane,fnamelabel{i})
    csvwrite(fname,recm(:,:,i,vindm(i)));
    figure(80+i)
    imshow(recm(:,:,i,vindm(i)),[])
    title(fname)
end

figure(100)
imshow(recm(:,:,5,vindm(5))-recm(:,:,6,vindm(6)),[])


%%

plane=720;
nshift=11; % number of tests
minshift=0; % in pixels
maxshift=2; % how far shall we shift
extrashift=linspace(minshift,maxshift,nshift);

sino=squeeze(d(:,plane,:));

v=nan(nshift,2);
baseshift=0*squeeze(T.Cen.fitshift(cas,rep,:));

rec=zeros(recsize,recsize,nshift);
MaskThresh=T.Cen.MaskThresh;

mask(:,:,1)=f.MakeTomoMask('par',f.fraccircshift(...
    squeeze(sino),-baseshift(plane)),recsize,ang,detPitch,src,det,...
    MaskThresh,0.000005,[150 150]);

mask(:,:,2)=f.MakeTomoMask('fan',f.fraccircshift(...
    squeeze(sino),-baseshift(plane)),recsize,ang,detPitch,src,det,...
    MaskThresh,0.000005,[150 150]);
f.CheckMMask(mask(:,:,1),10,10)
title('parallel')
pause(1)
f.CheckMMask(mask(:,:,2),10,10)
title('fan')


q=1;
tic

for i=1:nshift
    disp(i)
    sinogram=f.fraccircshift(sino,-(baseshift(plane)+extrashift(i)))';
    rec(:,:,i,1)=a.FBPexpl(sinogram,recsize,ang);
    rec(:,:,i,2)=a.FBPexplFan(sinogram,...
        recsize,ang,detPitch,src,det);
end

% plot
showim=[];
showimt=[];
for i=1:nshift
        rt1=rec(:,:,i,1);
    rt1m=mean(rt1(:));
        rt2=rec(:,:,i,2);
    rt2m=mean(rt2(:));
    for j=1:2
        [v(i,j),~]=qc.VarianceQuality(rec(:,:,i,j),mask(:,:,j));
    end
end
[vmax,vind]=max(v);
%
for i=1:nshift
    
    showimt=cat(1,rec(20:90,60:130,i,1)/rt1m,rec(20:90,60:130,i,2)/rt2m);
    if i==vind(2)
       showimt(1:10,1:10)=max(showimt(:)) ;
    end
        
    showim=cat(2,showim,showimt);
end
%

figure(1);clf;
imshow(showim,[]);
title('center-shift scan for parallel (top) and fan (bottom). white square marks best one')
figure(2);clf;
plot(extrashift,v(:,1)/max(v(:,1)),'rx','displayname','Par')
hold on
plot(extrashift,v(:,2)/max(v(:,2)),'+b','displayname','fan')
legend('Location','best')
xlabel('center-shift')
ylabel('quality criterion, normalized')
 title('center shift scan')
grid on
%
figure(3);clf % tomogram with path
imshow(rec(:,:,vind(2),2),[])
hold on
profx=[31,54];
profy=[86,79];
plot(profx,profy,'r')
title('best fan tomo with path')

figure(5);clf % tomogram with path
imshow(rec(:,:,vind(1),1),[])
hold on
profxp=[26,50];
profyp=[84,77];
% profxp=profx;
% profy=profy;

plot(profxp,profyp,'g')
title('best parallel tomo with path')

% path profiles
figure(7);clf

ax1=gca;
title('fan')
hold on
grid on

col=hsv(nshift);
for i=1:nshift
    
    p=plot(improfile(rec(:,:,i,2),profx,profy),'color',col(i,:),'DisplayName',sprintf('%.2f',extrashift(i)));
    if i==vind
        p.LineWidth=2;
        p.Color=[0 0 0];
        p.DisplayName=strcat(p.DisplayName,' (best)');
    end
end
legend('Location','best')
xlabel('path position')
ylabel('quality criterion')
figure(8);clf
title('parallel')
hold on
grid on

for i=1:nshift    
    p=plot(improfile(rec(:,:,i,1),profxp,profyp),'color',col(i,:),'DisplayName',sprintf('%.2f',extrashift(i)));
    if i==vind
        p.LineWidth=3;
        p.Color=[0 0 0];
        p.DisplayName=strcat(p.DisplayName,' (best)');
    end

end
xlabel('path position')
ylabel('quality criterion')
legend('Location','best')

figure(9);clf
title('fan - parallel comparison, nomalized')
hold on
grid on
p1=improfile(rec(:,:,vind(2),2),profx,profy);
p2=improfile(rec(:,:,vind(1),1),profxp,profyp);

plot(p1/max(p1),'r','displayname','best fan')
plot(p2/max(p2),'b','displayname','best parallel')
xlabel('path position')

legend('Location','best')

























