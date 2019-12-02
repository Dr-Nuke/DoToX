function [ output_args ] = f42_ArcLineIntercept( l1,l2,c,r,p0,dp)
% this function checks if a given Arc has antersections with a given line
% l1= [x,y] of line start point
% l2= [x,y] of line end point
% c = [x,y] of arc center point
% r = radius of arc
% p0= start angle of arc
% dp= angular size of the arc

dl=l2.-l1; % 

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

