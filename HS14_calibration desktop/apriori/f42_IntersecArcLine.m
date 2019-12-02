function [ output_args ] = f42_IntersecArcLine( P,i,x1,y1,x2,y2 )
% finds the intersection of an Arc an a line

% first make compute the linecirc arguments
a=(y2-y1)/(x2-x1);  %slope
b=y1-a*x1;          %y-axis intercept

% run the linecirc function
[x,y] = linecirc(a,b,P.c(i,1),P.c(i,2),P.r(i))

%then check if it is within the angular range


end

