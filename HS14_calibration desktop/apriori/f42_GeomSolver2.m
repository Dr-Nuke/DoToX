function [p3 x4 y4 t] = f42_GeomSolver2()
%this function computes the channel geometry values for the 
% two points in the D2O channel
% required for the analytical channel model
%% parameters
x3=6.7;
y3=6.7;
r3=4.14;

r4=2;
p4=pi/4;

x6=6.7+1/sqrt(2);
y6=6.7+1/sqrt(2);

%% Equation system (derived from intersecting the arcs + tangent requirement)
syms p3 x4 y4 t;
eqns = [x3+r3*cos(p3) == x4+r4*cos(p3),...
        y3+r3*sin(p3) == y4+r4*sin(p3),...
        x6+t*1/sqrt(2)== x4+r4*cos(p4),...
        y6-t*1/sqrt(2)== y4+r4*sin(p4)];

% numerical solver
T=vpasolve(eqns, [p3 x4 y4 t], [5.0115 7.331 4.655 1.8920]);

p3=double(T.p3);
x4=double(T.x4);
y4=double(T.y4);
t=double(T.t);

end


% cheat sheet
% r=[13.604,2,5.14,2,13.604];                             % radii
% ai=[7.07,125.34,-174.83,125.34,7.07];ai=ai/360*(2*pi);  % angle increment in rad
% a0=[0,7.07,-47.59,-42.41,90-7.07];a0=a0/360*(2*pi);     % statr angles in rad
% cxy=[0,11.516,6.7,1.429,0;...                           % center points
%     0,1.429,6.7,11.516,0];