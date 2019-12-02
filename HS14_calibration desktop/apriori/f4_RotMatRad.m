function [ R ] = f4_RotMatRad( a )
% creates the 2d rotation matrix from an angle

R=[cos(a),-sin(a);...
   sin(a), cos(a)];

end

