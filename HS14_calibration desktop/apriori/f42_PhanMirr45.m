function [P] = f42_PhanMirr45(P)
% this function mirrors the phantom geometry at the 45° line
% is specific for the creating of the channel hantom. does not handle 
% special cases of slopes being zero & inf 

T=[0,1;1,0]; %transform matrix

%     P.n=[];            % segment number
%     P.CL=[];            % Circle- Line indicator: equals 1 if the segment is a circle   
    P.c=(T*P.c')';       % circle cter
%    P.r=[];           % radius
    P.p0=pi/2-P.p0;          % circle segment start angle
    P.dp=-P.dp ;        % circle segment angular size
%    P.reg=[];
    P.Lx1=(T*P.Lx1')'; %starting points
    P.Lx2=(T*P.Lx2')'; % ending points

    P.dx=(T*P.dx')';
    P.a=1./P.a;
    P.b=(P.Lx1(:,2)-(P.a').*P.Lx1(:,1))' ; %b=y-ax

end

