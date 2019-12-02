% check out the "transfer fuction" of forward and backward projections

rs=269; % canvas size / recon size
phan=zeros(rs);

%make phantom
phan(100:120,150:170)=1;
phan(150:170,110:130)=0.5;



%check phantom
figure(373);clf
imshow(phan',[]);set(gca,'YDir','normal')
axis on
title(sprintf('phantom, imsum=%.2f',sum(phan(:))))

%make angle list
angl=linspace(0,360,1357+1);
angl(end)=[];

%maks sinogram
sino = radon(phan,angl);

%check sino
figure(374);clf
imshow(sino,[]);set(gca,'YDir','normal')
colorbar()
axis on
title('sinogram')

% reconstruc
recon=iradon(sino,angl,'spline','Hann',1,rs);

recon2=a.FBPexplFan(...
                        sino',...
                        rs,deg2rad(angl),0.127,1000000,1);

%check the recon
figure(375);clf
imshow(recon',[]);set(gca,'YDir','normal')
axis on
colorbar
title(sprintf('recon irad matlab, imsum=%.2f',sum(recon(:))))

figure(376);clf
imshow(recon2',[]);set(gca,'YDir','normal')
axis on
colorbar
title(sprintf('recon astra fan GPU, imsum=%.2f',sum(recon(:))))

figure(377);clf;
histogram(recon2)