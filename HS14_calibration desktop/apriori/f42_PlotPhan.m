function f42_PlotPhan(P,ax,col,ind_d)
% makes a plot of the Ptom with correct segment heights
% f42_PlotPhan(P,ax,col,ind_d)
% P = the Phantom 
% c = color specifier
% ax = axis handle
% col = rgb color vector
% ind_d = draw indicator: vector specifying the content to be plotted.
%       1 = lines, 2= start points, 3 = end points, 4 = centers



ax=gca;

hold on


for i=1:length(P.n)
    if P.CL(i)==1 %if the secment is an arc
        f4_PlotArc(P.c(i,:),P.r(i),P.p0(i),P.dp(i),360,ax,col,ind_d)

        
    else % else if it is a straight line
        f4_PlotLine(P.Lx1(i,:),P.Lx2(i,:),ax,col,ind_d)
        
    end
end
        
        



grid on
xlabel('x')
ylabel('y')

axis equal


end
