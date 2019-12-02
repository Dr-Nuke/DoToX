close all
a=squeeze(d2o_corr(:,50,:));
imbo3(a,1);
b=mean(a,2);
figure(2);clf; plot(b,'.'); grid on; hold on

%%
c= smooth(b,0.1,'loess')
plot(c,'k')

%%
figure(3)
plot(b./c,'.')

%%
corr=repmat(c./b,1,375);
aa=a.*corr;
imbo3(aa,4);