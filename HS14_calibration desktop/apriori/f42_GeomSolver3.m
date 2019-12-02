function [x5 y5 p5 p45] = f42_GeomSolver3(x2,y2,x4,y4)
%this function computes the channel geometry values for the 
% two points in the D2O channel
% required for the analytical channel model
%% parameters

r23=3;
r45=3;
r5=2;



%% Equation system (derived from intersecting the arcs + tangent requirement)
syms x5 y5 p5 p45;
eqns = [x2+r23*cos(p5-pi)   == x5+r5*cos(p5),...
        y2+r23*sin(p5-pi)   == y5+r5*sin(p5),...
        x4+r45*cos(p45)     == x5+r5*cos(p45+pi),...
        y4+r45*sin(p45)     == y5+r5*sin(p45+pi)];

% numerical solver
T=vpasolve(eqns, [x5 y5 p5 p45], [12.015 6.404 4.6124 0.3573]);

x5=double(T.x5);
y5=double(T.y5);
p5=double(T.p5);
p45=double(T.p45);

end


% cheat sheet
% r=[13.604,2,5.14,2,13.604];                             % radii
% ai=[7.07,125.34,-174.83,125.34,7.07];ai=ai/360*(2*pi);  % angle increment in rad
% a0=[0,7.07,-47.59,-42.41,90-7.07];a0=a0/360*(2*pi);     % statr angles in rad
% cxy=[0,11.516,6.7,1.429,0;...                           % center points
%     0,1.429,6.7,11.516,0];