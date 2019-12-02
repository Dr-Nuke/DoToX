function [P] = f42_PhanMove(P,x,y)
% moves the Phantom according ro the input in the 2d plane
% x = x offset
% y = y offset

P.c(:,1)=P.c(:,1)+x;
P.c(:,2)=P.c(:,2)+y;

P.Lx1(:,1)=P.Lx1(:,1)+x;
P.Lx1(:,2)=P.Lx1(:,2)+y;

P.Lx2(:,1)=P.Lx2(:,1)+x;
P.Lx2(:,2)=P.Lx2(:,2)+y;

P.b=P.b+y-P.a*x;

end

