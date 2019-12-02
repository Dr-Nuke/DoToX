function [P] = f42_PhanZoom(P,z)
% zooms the phantom by a factor z
% P = the phantom
% z = zoom facotr

P.c=P.c*z;
P.r=P.r*z;
P.Lx1=P.Lx1*z;
P.Lx2=P.Lx2*z;
P.dx=P.dx*z;
P.b=P.b*z;

end