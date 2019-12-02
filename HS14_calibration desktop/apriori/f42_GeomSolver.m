function [ph1,ph2,x3,y3] = f42_GeomSolver()
%this function computes the channel geometry values for the 
% center point of the r=2 arcs and the according angles,
% both required for the analytical channel model
%% parameters
x1=0;
y1=0;
r1=13.604;

x2=6.7;
y2=6.7;
r2=5.14;

r3=2;

%% Equation system (derived from intersecting the arcs + tangent requirement)
syms p1 p2 x3 y3;
eqns = [x1+r1*cos(p1) == x3+r3*cos(p1),...
        y1+r1*sin(p1) == y3+r3*sin(p1),...
        x2+r2*cos(p2) == x3+r3*cos(p2+pi),...
        y2+r2*sin(p2) == y3+r3*sin(p2+pi)];

% numerical solver
T=vpasolve(eqns, [p1 p2 x3 y3], [0.1234 5.4526 11.516 1.429]);

ph1=double(T.p1);
ph2=double(T.p2);
x3=double(T.x3);
y3=double(T.y3);

end


% cheat sheet
% r=[13.604,2,5.14,2,13.604];                             % radii
% ai=[7.07,125.34,-174.83,125.34,7.07];ai=ai/360*(2*pi);  % angle increment in rad
% a0=[0,7.07,-47.59,-42.41,90-7.07];a0=a0/360*(2*pi);     % statr angles in rad
% cxy=[0,11.516,6.7,1.429,0;...                           % center points
%     0,1.429,6.7,11.516,0];