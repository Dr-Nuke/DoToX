ind=[1,6,7,8]
x=[474,605];
y=[585,1061];
yplane=200
prof(1,:)=improfile(squeeze(d1(:,yplane,:)),x,y);
prof(6,:)=improfile(squeeze(d6(:,yplane,:)),x,y);
prof(7,:)=improfile(squeeze(d7(:,yplane,:)),x,y);
prof(8,:)=improfile(squeeze(d8(:,yplane,:)),x,y);

fig=figure(12);clf;
ax=gca();
hold on;
for i=ind
    p(i)=plot(prof(i,:),'displayname',num2str(i));
    
end
legend()

dm(1)=mean(mean(mean(d1(470:end,130:1000,:),3),2));
dm(6)=mean(mean(mean(d6(470:end,130:1000,:),3),2));
dm(7)=mean(mean(mean(d7(470:end,130:1000,:),3),2));
dm(8)=mean(mean(mean(d8(470:end,130:1000,:),3),2));

figure(13);clf;
bar(dm)

%%

dd(6)=mean(mean(mean(div6(470:end,130:1000,:),3),2));
dd(7)=mean(mean(mean(div7(470:end,130:1000,:),3),2));
dd(8)=mean(mean(mean(div8(470:end,130:1000,:),3),2));

figure(14);clf;
bar(dd)