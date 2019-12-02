figure(2)
clf
plot(phan.c,'LineWidth',2)
hold on
plot(xtot(1,:),'.','Markersize',10)
plot(xtot(2,:),'x','Markersize',10)

legend({'phantom','zero noise','30% noise'},'Location','southeast')

xlabel('segment number')
ylabel('attenuation')
set(gca,'Fontsize',15)