%%

load('/media/data/cbolesch/Tomo_HS14/processed/block_4add_cl3_.mat')
load('/media/data/cbolesch/Tomo_HS14/processed/block_4add_d2o_.mat')



cl3=block_4add_cl3_;
clear block_4add_cl3_
d2o=block_4add_d2o_;
clear block_4add_d2o_
%%


k=350;
sinocl3=squeeze(cl3(:,k,:));
sinod2o=squeeze(d2o(:,k,:));

figure(10); clf
ax1=subplot(3,1,1);
imbo4(sinocl3);

ax2=subplot(3,1,2);
imbo4(sinod2o);

ax3=subplot(3,1,3);
loga=-log(sinocl3./sinod2o);
imbo4(loga);


linkaxes([ax1,ax2,ax3],'xy')