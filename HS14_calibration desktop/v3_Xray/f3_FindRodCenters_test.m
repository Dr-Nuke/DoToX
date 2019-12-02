% tests the center point fit function. assumes F already exists

figure(23)
clf
range=100:420;

F=f3_FindRodCenters(F,range)



for i = 1:4
    subplot(2,2,i)
    for j=1:2
        plot(F.c(:,i,j),'b')
        hold on
        plot(F.cfit(:,i,j),'r--')
    end
end

        %%
figure(25)
clf
gca
hold on
for i=1:4
    
        %plot3((F.cfit(:,i,1)-0.99*F.cfit(1,i,1))/res,(F.cfit(:,i,2)-0.99*F.cfit(1,i,2))/res,1:490)
        %plot3((F.cfit(1,i,1)-0.99*F.cfit(1,i,1))/res,(F.cfit(1,i,2)-0.99*F.cfit(1,i,2))/res,1:490,'o')
        plot3((F.cfit(:,i,1))/res,(F.cfit(:,i,2))/res,1:490)
        plot3((F.cfit(1,i,1))/res,(F.cfit(1,i,2))/res,1:490,'o')
    
end

xlabel('x')
ylabel('y')
zlabel('z')
grid on
