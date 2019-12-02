i=2;
j=1;

figure(1)
hold on
plot([A.cs.center(i,1),A.paths.endpoints{i}(j,1)],...
            [A.cs.center(i,2),A.paths.endpoints{i}(j,2)],'or')

dist=sqrt((cx-A.cs.center(4,1)).^2+(cy-A.cs.center(4,2)).^2)/res;

figure(2)
hold on
col =hsv(10)
for j=1:10
    
plot(dist,squeeze(A.profiles.profile(i,j,:)),'Color',col(j,:))
end
line([5.14,5.14],[0,1], 'Color',[0 0 0])
line([4.14,4.14],[0,1], 'Color',[0 0 0])