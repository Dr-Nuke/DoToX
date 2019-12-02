function [P] = f42_RotPhan(P,a)
% rotates the phantom in a given angle
% P  = original phantom
% a     = rotation angle


R=f4_RotMatRad(a);         % create 2d rotation matrix
P.c=(R*P.c')';              % Rotate arc center  point coordinates
P.p0=mod(P.p0+a,2*pi);  % shift starting angles
P.Lx1=(R*P.Lx1')';            % rotate the start and end points of the line segments
P.Lx2=(R*P.Lx2')';

P.dx=(R*P.dx')';

P.a=tan((atan(P.a))+a);

P.b=(P.Lx1(:,2)-(P.a').*P.Lx1(:,1))' ; %b=y-ax
end

