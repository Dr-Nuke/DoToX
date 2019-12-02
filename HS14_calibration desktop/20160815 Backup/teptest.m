figure(5);cla; hold on
figure(6);cla; hold on
h=[550:100:1950]
hh=length(h);
c=hsv(hh);



for i=1:hh
    figure(5)
    plot(G.p(h(i)).pin(4).thickness,'color',c(i,:))
    figure(6)
    plot(i,G.p(h(i)).pin(1).th_mean,'ro')
    plot(i,G.p(h(i)).pin(2).th_mean,'bo')
    plot(i,G.p(h(i)).pin(3).th_mean,'go')
    plot(i,G.p(h(i)).pin(4).th_mean,'mo')
end
    %%
figure(7)
cla
hold on
num=length(G.p(950).pin(1).prof(:,1));
c=hsv(num);
for i=1:5:num
    plot(dist,G.p(950).pin(1).prof(i,:),'Color',c(i,:));
end

%%

figure(8)
hold on
c=hsv(hh)
for i=1:hh
    figure(5)
    plot(G.p(h(i)).pin(1).thickness,'color',c(i,:))
end