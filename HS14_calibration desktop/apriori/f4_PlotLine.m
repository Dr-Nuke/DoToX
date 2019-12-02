function  f4_PlotLine(x1,x2,ax,col,ind_d)
% plots a line with special end points
% ind_d = draw indicator: vector specifying the content to be plotted.
%       1 = lines, 2= start points, 3 = end points, 4 = centers

if any(1==ind_d)
line(ax,[x1(1),x2(1)],[x1(2),x2(2)],'Color',col)
end
if any(2==ind_d)
plot(ax,x1(1),x1(2),'o','Color',col)
end
if any(3==ind_d)
plot(ax,x2(1),x2(2),'x','Color',col)
end

% xmax=max(x(:));
% xmin=min(x(:));
% ymax=max(y(:));
% ymin=min(y(:));
% grid on
% axis([xmin xmax ymin ymax])
% axis equal


end

