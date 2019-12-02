figure(1)


for i =1:4
subplot(4,1,i)
imbo4(LFT(:,i,:))
ylabel(sprintf('angle [deg]'))
if i==1,
    title(sprintf('non-scaled LFT \n rod %d',i))
else
    title(sprintf('rod %d',i))
end
set(gca,'YTick', [1,178]);
set(gca,'YTickLabel', [0,90]);
%ax.YTickLabel = {'0','90'};

end
xlabel('rod height [pixel]')


figure(2)
hold on
c=hsv(4)
for i =1:4
    plot(mean(squeeze(LFT(:,i,:)),1),'Color',c(i,:))
end
legend()
    
figure(3)
hold on
c=hsv(4)
for i =1:4
    plot(mean(squeeze(LFT(:,i,:)),2),'Color',c(i,:))
end
legend()
    