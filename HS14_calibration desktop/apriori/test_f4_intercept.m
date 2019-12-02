% test this function.

% finds the interception of two lines, defined by their end point coordinates 
% xy10 = [x,y] of the 1st point of the 1st line
% xy11 = [x,y] of the 2nd point of the 1st line
% xy20 = [x,y] of the 1st point of the 2nd line
% xy21 = [x,y] of the 2nd point of the 2nd line


xy10=[0,0];
xy11=[1,0];
xy20=[0.5,-0.5];
xy21=[0.5,0.5];

[t,cosphi]=f4_intercept(xy10,xy11,xy20,xy21);

int1=xy10+t(1)*(xy11-xy10);
int2=xy20+t(2)*(xy21-xy20);


figure(1)
clf
plot(xy10(1),xy10(2),'x','DisplayName','x10')
hold on
plot(xy11(1),xy11(2),'x','DisplayName','x11')
plot([xy10(1),xy11(1)],[xy10(2),xy11(2)],'DisplayName','line 1')

plot(xy20(1),xy20(2),'+','DisplayName','x20')
plot(xy21(1),xy21(2),'+','DisplayName','x21')
plot([xy20(1),xy21(1)],[xy20(2),xy21(2)],'DisplayName','line 2')
plot(int1(1),int1(2),'o','DisplayName','intersect 1')
plot(int2(1),int2(2),'o','DisplayName','intersect 2','Markersize',10)
legend('show')
disp([cosphi])

xlim([-3,3])
ylim([-3,3])

grid on
title(sprintf('intersection at t1 = %.2f and 2t = %.2f \n cosphi = %f',t(1),t(2),cosphi))