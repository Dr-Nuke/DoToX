function [phan] = f4_PhanRotate(phan,a)
% rotates the phantom in a given angle
% phan  = original phantom
% a     = rotation angle


R=f4_RotMat(a);
phan.xy=R*phan.xy;

end

