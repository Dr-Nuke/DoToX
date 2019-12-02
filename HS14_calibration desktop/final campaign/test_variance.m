cas=1;
rep=1;

if 0
    fname=sprintf('%s%02d_%02d.mat','1_corr_',T.cas,T.rep);
    fpath=sprintf('%s%s',T.d.DataPath,fname);
    d=f.BoLoad(fpath,T);
    
    
    % image correction: pixels & filter
    corr=f.ImCorrFilt(mean(d,3),T);
    % apply image correction
    d=f.ApplyImCorr(d,corr,T);
    
    % crop x,z
    d=d(T.Cen.range,:,T.q360.startframe:(...
        T.q360.startframe+T.q360.nFrames(1,1)-1));
    
    % quadrant-rotate, fliplr, log
    d=flip(-log(circshift(d,T.Cen.quadrot,3)),3);
    
    % find the rough centering for each frame ~+-2pix
    [T.Cen.fit{cas,rep},T.Cen.fitshift(cas,rep,:),T.Cen.centshift(cas,rep,:)]=...
        f.FindCentering(T,d);
    
    
    T.fnames.corr2{cas,rep}=f.BoSave3(d,'1_cor2_',T);
else
    d=f.BoLoad(T.fnames.corr2{1,1},T);
end
%%
recsize=size(d,1);
sinsh=T.Rec.SinoRotShift;
ang=deg2rad(T.Rec.angles);
recsize=T.Rec.recsize;
detPitch=T.Rec.detPitch;
src=T.Rec.src;
det=T.Rec.det;
niter=T.Cen.niter;
BS=T.Raw.BS;
xmin=T.Cen.xmin;
xmax=T.Cen.xmax;
MaskThresh=T.Cen.MaskThresh;

ang=deg2rad(T.Rec.angles);
baseshift=squeeze(T.Cen.fitshift(cas,rep,:));
MaskThresh=T.Cen.MaskThresh;

T.Cen.MMask=f.MakeTomoMask('par',f.fraccircshift(...
    squeeze(d(:,T.Cen.MaskPlane,:)),...
    -baseshift(T.Cen.MaskPlane)),recsize,ang,detPitch,src,det,MaskThresh,...
    0.000005,[150 150]);
f.CheckMMask(T.Cen.MMask,cas,rep)

MMask=T.Cen.MMask;

niter=20;
v=nan(niter,BS(2));
xlist=zeros(niter,BS(2));
eshift=zeros(1,BS(2));
recon=zeros(recsize,recsize,BS(2),'single');

fprintf('planes: ')
for plane=720% 10:1017%plane=[50,1000]%
    
    f.f_BoCount(plane,100,10,5)
    
    % check out sinogram
    sino=squeeze(d(:,plane,:));
    
    % center exact
    [eshift(plane),xlist(:,plane),...
        v(:,plane)]=f.fibo(sino,0,...
        xmin,xmax,niter,0.007,1,ang,recsize,detPitch,src,det,MMask);
    
    
    
    % reconstruct
    recon(:,:,plane)=a.FBPexplFan(...
        f.fraccircshift(sino,-(eshift(plane)+baseshift(plane)))',...
        recsize,ang,detPitch,src,det);
    
    
    f.CheckRecon(recon(:,:,plane)...
        ,v(:,plane),xlist(:,plane)...
        ,cas,rep,plane,xmin,xmax)

end
