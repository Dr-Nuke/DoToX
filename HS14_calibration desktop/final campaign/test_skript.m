

plane=1000;

crange=[-1,9]*0.001
fig(1)=figure(1);clf
set(fig(1),'Position',[81 574 493 412])
ax(1)=imshow(squeeze(recona(:,:,plane,1)),[])
title('empty absolute')

fig(2)=figure(2);clf
set(fig(2),'Position',[81 574 493 412])
ax(1)=imshow(squeeze(recona(:,:,plane,1)),[])
title('empty absolute')
