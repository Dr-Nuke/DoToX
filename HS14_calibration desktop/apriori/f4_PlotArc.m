function  f4_PlotArc(c,r,ph0,dph,n,ax,col,ind_d)
% plots an Arc, end and center points as specified
% c = center point [x,y]
% r = radius
% ph0 = start angle (phi_zero)
% dph = delta phi, the arc angle
% n = the point density, => n/360 points per full circle
% ax = axis handle
% col = [r,g,b] vector for coloring
% ind_d = draw indicator: vector specifying the content to be plotted.
%       1 = lines, 2= start points, 3 = end points, 4 = centers

a=linspace(ph0,ph0+dph,round(n*abs(dph)/(2*pi)));

x=r*cos(a)+c(1);
y=r*sin(a)+c(2);

if any(1==ind_d)
    plot(ax,x,y,'Color',col)
end
if any(2==ind_d)
plot(ax,x([1]),y([1]),'o','Color',col)
end
if any(3==ind_d)
plot(ax,x([end]),y([end]),'x','Color',col)
end
if any(4==ind_d)
plot(ax,c(1),c(2),'x','Color',col)
end

% xmax=max(x(:));
% xmin=min(x(:));
% ymax=max(y(:));
% ymin=min(y(:));
% grid on
% axis([xmin xmax ymin ymax])
% axis equal


end

