im=squeeze(mean(d,3));
med=medfilt2(im,[3 3]);
med=imgaussfilt(im, 1);

diff=abs(med-im);
figure(400);clf;
nplot=4;
crange=[0.55,0.65];

ax(1)=subplot(2,2,1);
imshow(im,crange)
xlim()

ax(2)=subplot(2,2,2);
imshow(med,crange)

ax(3)=subplot(2,2,3);
imshow(diff,[])
ax(4)=subplot(2,2,4);

imshow(abs(im-med)>0.01)
linkaxes(ax)

%%
figure(401);clf
a=mean(d,3);
imshow(a,[0.5,0.6])

xlim([465,485]);
ylim([285,305]);

figure(402);clf
ax(5)=subplot(2,1,1);
hold on
grid on


xstart=292;
ystart=470;
for i=1:5
    ind=xstart+i-1
    p(i)=plot(a(ind,:),'-x','displayname',num2str(ind));
end
xlim([465,485])
legend(p)


ax(6)=subplot(2,1,2);
hold on
grid on
for i=1:8
    ind=ystart+i-1
    pp(i)=plot(a(:,ind),'-x','displayname',num2str(ind));
end
xlim([285,305])
legend()

%% => enter manual data here
xcor=293:295;
ycor=472:475;

xbase=290:298;
ybase=468:480;

mask=zeros(size(a));
mask(xbase,ybase)=1;
mask(xcor,ycor)=0;
figure(402);clf;
axx(1)=subplot(2,1,1);

c=[0.55,0.65];
imshow(a,c);
axx(2)=subplot(2,1,2);
imshow(imgaussfilt(a, 1),c)
linkaxes(axx)


