function [ R ] = f4_RotMat( a )
% creates the 2d rotation matrix from an angle
a=deg2rad(a);
R=[cos(a),-sin(a);...
   sin(a), cos(a)];

end

