figure(100)
cla
hold on
c=hsv(hh);
for i=1:hh
    plot(cumsum(mean(G.p(h(i)).pin(3).prof(:,1:end))),'Color',c(i,:))
end;;