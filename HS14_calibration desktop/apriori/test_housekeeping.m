figure(3);
clf
imshow(im); % image of the door
hold on;
x=[789,796];    % 
y=[894,182];
line(x,y,'LineWidth',2);
y2=[894,37];
line(x,y2,'Color','red')
d(1)=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
d(2)=sqrt((x(2)-x(1))^2+(y2(2)-y2(1))^2);

d=d/d(1)*2