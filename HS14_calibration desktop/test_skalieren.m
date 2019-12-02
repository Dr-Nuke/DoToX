figure(25)
cla
plot(mean(squeeze(LFT(:,1,:)),2))
hold on
scaling=0.0002
for i=0:800
    plot(i,f_plane2lft(i+550)*scaling,'r')
end
    