% assume we have planes 50 100 500 and 1000 in both recon and recon2
% (classic and absolute-diff recon)

% 1 make raw figures

reconsize=size(recon,1);

fhi=412;
fwi=493;
fx=[4 0 1 2 ]*fwi;

crange=[-0.003 0.005;
-1.5e-04 0.001;
-1.5e-04 0.001;
-1.5e-04 0.001];


plane=[50 100 500 1000];
pla=2;
coords=[32,279;
    79,234];

for cas = 1:4
    
fh(cas)=figure(cas);clf;
fh(cas).Position=[fx(cas) 1 fwi fhi];
ax(cas)=imshow(recon(:,:,plane(pla),cas)',[crange(cas,:)]);
set(gca,'YDir','normal');
hold on
plot(coords(1,:),coords(2,:),'color','r');
title(sprintf('A: -ln(F/nF)[%d,%d]',cas,plane(pla)))
end

for cas = 1:4
    
fh(cas+4)=figure(cas+4);clf;
fh(cas+4).Position=[fx(cas) 550 fwi fhi];
ax(cas)=imshow((recon2(:,:,plane(pla),cas)-...
    recon2(:,:,plane(pla),1))',[crange(cas,:)]);
set(gca,'YDir','normal');
hold on
line(coords(1,:),coords(2,:),'color','r');
title(sprintf('B: -ln(F)-ln(nF) [%d,%d]',cas,plane(pla)))
end

%%

figure(10);clf

prof1=improfile(squeeze(recon(:,:,plane(pla),4)),coords(1,:),coords(2,:));
prof2=improfile(squeeze(recon2(:,:,plane(pla),4))...
    -recon2(:,:,plane(pla),1),coords(1,:),coords(2,:));

plot(prof1,'Displayname','A')
hold on
plot(prof2,'Displayname','b')
grid on
legend