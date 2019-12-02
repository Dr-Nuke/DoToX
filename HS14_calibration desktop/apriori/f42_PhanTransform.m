function [P] = f42_PhanTransform(P,T)
%transorms the coordinates of a Phantom according to a transform Matrix
% P = phantom
% T = transform Matrix



%     P.n=[];            % segment number
%     P.CL=[];            % Circle- Line indicator: equals 1 if the segment is a circle   
    P.c=(T*P.c')';       % circle cter
%    P.r=[];           % radius
    P.p0=[];          % circle segment start angle
    P.dp=[] ;        % circle segment angular size
    P.reg=[];
    P.Lx1=[]; %starting points
    P.Lx2=[]; % ending points


end

