function [x,y] = f_LinGen(center, a_start,da,n,r)
%F_LINGEN generates points that have acertain distance and angle towards a
% common central point. 
% center    = [x,y] the central point of reference
% a_start   = ang a start angle in degree
% da        = ang2 the angle range covered. negative for clockwise, positive for
% counter clockwise
% n         = n number of line end points generated, equally distributed over [a_start,
% a_start+da]
% r         = radius for the generated points

%construct angles
alpha= degtorad(a_start+linspace(0,da,n));


x=center(1)+cos(alpha)*r;
y=center(2)+sin(alpha)*r;

end

