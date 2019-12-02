%% plotscriiiiiipt!

%% Preprocessing
% delete(gcp('nocreate'));
% parpool(7);
% preprocess
for cas=1%1:T.d.ncas %iterate cases
    if any(T.d.nframes(cas,:)~=0) % skip if no files available
        T.cas=cas; % copy iterater to T struct
        for rep=1% iterate repetitions
            d=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
            if T.d.nframes(cas,rep)~=0 % skip if no measurement available
                T.rep=rep; % copy iterater to T struct
                fprintf('%s preprocessing...\n',f.CommLineStart(cas,rep)); %debug
                if 0 % skip some data processing steps if we did them already
                    % 1 read in
                    d=f.ReadXrayTomo(T);
                    T.Raw.rad(:,:,cas,rep)=squeeze(d(:,:,1));
                    
                    % 2 find rot startt & stop and crop
                    T.Raw.CropRange(cas,rep,:,:)=f.FindCrop(T,...
                        squeeze(d(:,T.Raw.sinoyplanes(1),:)));
                    
                    % make 360° check plot
                    T.q360.nFrames(cas,rep)=...
                        f.frames360(squeeze(d(:,T.Raw.sinoyplanes(1),:)),T);
                    
                    % plot the graph
                    f.CheckSino1(T,d,10);
                    
                    % crop & copy
                    [d,T.Raw.nframes(cas,rep)]=f.SinoCrop(d,T);
                    
                    % 3 dose correction, some plots
                    [d,T.dose.doses(cas,rep,:),dosemin]=f.DoseCorrect(d,T);
                    f.DoseMaskCheck(d(:,:,1),T,40);
                    f.DoseCorrectCheck1(T,dosemin,50);
                    f.DoseCorrectCheck2(d,T,70);
                    clear dosemin
                    
                    % 4 beam correct
                    
                    [d,T.fit.fits{cas,rep},T.fit.fitdose{cas,rep},im_old]=...
                        f.BeamCorrect(d,T);
                    f.BeamCorrectCheck(d,T,im_old);
                    f.DoseCorrectCheck2(d,T,71);
                    clear im_old
                    
                    % save intermediate result
                    T.fnames.corr{cas,rep}=f.BoSave3(d,'1_corr_',T);
                    %error('lets stop here')
                    
                else
                    fname=sprintf('%s%02d_%02d.mat','1_corr_',T.cas,T.rep);
                    fpath=sprintf('%s%s',T.d.DataPath,fname);
                    d=f.BoLoad(fpath,T);
                end
                
                
                % image correction: pixels & filter
                corr=f.ImCorrFilt(mean(d,3),T);
                % apply image correction
                d=f.ApplyImCorr(d,corr,T);
                
                % crop x,z
                d=d(T.Cen.range,:,T.q360.startframe:(...
                    T.q360.startframe+T.q360.nFrames(1,1)-1));
                
                % quadrant-rotate, fliplr, log
                d=flip(-log(circshift(d,T.Cen.quadrot,3)),3);
                
%                 % find the rough centering for each frame ~+-2pix
%                 [T.Cen.fit{cas,rep},T.Cen.fitshift(cas,rep,:),T.Cen.centshift(cas,rep,:)]=...
%                     f.FindCentering(T,d);
%                 
                baseshift=squeeze(T.Cen.fitshift(cas,rep,:));

                    T.Cen.MMask(:,:,cas,rep)=f.MakeTomoMask('par',f.fraccircshift(...
                        squeeze(d(:,T.Cen.MaskPlane,:)),...
                        -baseshift(T.Cen.MaskPlane)),recsize,ang,detPitch,src,det,...
                        MaskThresh,0.000005,[150 150]);
                    f.CheckMMask(T.Cen.MMask(:,:,cas,rep),cas,rep)
                    % end
                MMask=T.Cen.MMask(:,:,cas,rep);
                
                v=nan(niter,BS(2));
                xlist=nan(niter,BS(2));
                eshift=nan(1,BS(2));
                
                
                
                
            end   % isempty(T.fnames.add{cas}) % skip if no files available
        end
    end
end


%%

plane=512;
% check out sinogram
sino=squeeze(d(:,plane,:));

% center exact
[eshift(plane),xlist(:,plane),...
    v(:,plane)]=f.fibo(sino,0,...
    xmin,xmax,niter,0.007,1,ang,recsize,detPitch,src,det,MMask);
eshift(plane)=round(eshift(plane),2);
%baseshift(plane)

% reconstruct
rec1=a.FBPexplFan(...
    f.fraccircshift(sino,-eshift(plane))',...
    recsize,ang,detPitch,src,det);
rec2=a.FBPexpl(...
    f.fraccircshift(sino,-eshift(plane))',...
    recsize,ang);


f.CheckRecon(recon(:,:,plane)...
    ,v(:,plane),xlist(:,plane)...
    ,cas,rep,plane,xmin,xmax,0)
%%
fh=figure(134);clf;
pub.BFfigure(fh,T.F)
rm1=mean(rec1(:));
rm2=mean(rec2(:));

winx=(20:70);
winy=(70:120);

im1=cat(2,rec1'/rm1,rec2'/rm2);
im2=cat(2,imresize(rec1(winx,winy)',[recsize,recsize])/rm1,...
          imresize(rec2(winx,winy)',[recsize,recsize])/rm2);

im3=cat(1,im2,im1);      

imshow(im3,[])
hold on
set(gca,'YDir','normal')
rectangle('Position',[winx(1),recsize+winy(1),length(winx),length(winy)],...
    'edgecolor',[1 1 1])
rectangle('Position',[recsize+winx(1),recsize+winy(1),length(winx),length(winy)],...
    'edgecolor',[1 1 1])
plot([0,winx(1)],recsize+[0,winy(1)],'w')
plot([winx(end),recsize],recsize+[winy(2),0],'w')

plot(recsize+[0,winx(1)],recsize+[0,winy(1)],'w')
plot(recsize+[winx(end),recsize],recsize+[winy(2),0],'w')

axis on
ax=gca;
pub.BFaxis(ax,T.F)
xticks(recsize*[0.5,1.5])
xticklabels({'fan-beam','parallel-beam'})
yticks([])


set(ax,'TickLength',[0 0])
fh.PaperSize=max(fh.PaperSize,fh.Position(3:4));
fh.PaperPosition(3:4)=max(fh.PaperPosition(3:4),fh.Position(3:4));


print(sprintf('%sParVsFan.pdf',T.Fig.saveto),'-dpdf')
open(sprintf('%sParVsFan.pdf',T.Fig.saveto))



%% pixel correction plot
d=f.BoLoad(fpath,T);
               
fid=openfig('Fig 090 01 01 Image filter.fig');


im=mean(d(200:450,6:1019,:),3);
im2=fid.Children(1).Children.CData(6:1019,:);
close(fid)
%%
fid2=figure(1234);clf
pub.BFfigure(fid2,T.F)
imshow(cat(2,im',(im2-1)*50+0.5),[0,1])
hold on
set(gca,'YDir','normal')

p1 = [340 465;
    454 558;
    342 570;
    320 355];                         % First Point
p2 = [15 20;
    15 20;   
    15,20
    -23,-5];
p1=p1-p2;% Second Point
                        % Difference

quiver(p1(:,1),p1(:,2),p2(:,1),p2(:,2),0,'color','k')

axis on
ax=gca;
pub.BFaxis(ax,T.F)
xticks([120,370])
xticklabels({'time average','filter kernel'})
yticks([])

set(ax,'TickLength',[0 0])
fid2.PaperSize=max(fid2.PaperSize,fid2.Position(3:4));
fid2.PaperPosition(3:4)=max(fid2.PaperPosition(3:4),fid2.Position(3:4));

print(sprintf('%sImfilter.pdf',T.Fig.saveto),'-dpdf')
open(sprintf('%sImfilter.pdf',T.Fig.saveto))

%% sinogram plot
d=f.BoLoad(fpath,T);
%%
fh=figure(567);
clf;
pub.BFfigure(fh,T.F)
imshow(-log(squeeze(d(200:450,100,:))),[])
set(gca,'YDir','normal')
% axis on
% xlh=xlabel('frames');
% pub.BFxlab(xlh,T.Fig)
% ylh=ylabel('detector pixels');
% pub.BFylab(ylh,T.Fig)
sx=16;
sy=3.5;
fh.Position(3:4)=[sx sy];
fh.PaperPosition=[0 0 sx sy];
fh.PaperSize=[sx sy];
print(sprintf('%sSino.pdf',T.Fig.saveto),'-dpdf')
open(sprintf('%sSino.pdf',T.Fig.saveto))
















