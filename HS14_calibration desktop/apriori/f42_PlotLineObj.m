function [ output_args ] = f42_PlotLineObj(l,ax,ind_d,col)

%plots a line object
% l = the line object
% ax = axis handle
% col = [r,g,b] vector for coloring
% ind_d = draw indicator: vector specifying the content to be plotted.
%       1 = lines, 2= start points, 3 = end points, 4 = centers
hold on

for i = 1:length(l.a)

    if any(1==ind_d)
        line(ax,[l.xy1(i,1),l.xy2(i,1)],[l.xy1(i,2),l.xy2(i,2)],'Color',col)
    end
    if any(2==ind_d)
       plot(ax,l.xy1(i,1),l.xy1(i,2),'o','Color',col)
    end
    if any(3==ind_d)
       plot(ax,l.xy2(i,1),l.xy2(i,2),'x','Color',col)
    end



end

