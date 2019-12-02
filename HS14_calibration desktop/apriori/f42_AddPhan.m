function [P] = f42_AddPhan(varargin)
%% ads or concatenates two phantoms


    P.n=[];            % segment number
    P.CL=[];            % Circle- Line indicator: equals 1 if the segment is a circle   
    P.c=[];       % circle cter
    P.r=[];           % radius
    P.p0=[];          % circle segment start angle
    P.dp=[] ;        % circle segment angular size
    P.reg=[];
    P.Lx1=[]; %starting points
    P.Lx2=[]; % ending points
    
    P.dx=[]; % differential vector
    P.a=[];
    P.b=[];

for i =1:nargin
    P.n=horzcat(P.n,varargin{i}.n);
    P.CL=horzcat(P.CL,varargin{i}.CL);
    P.c=vertcat(P.c,varargin{i}.c);
    P.r=horzcat(P.r,varargin{i}.r);
    P.p0=horzcat(P.p0,varargin{i}.p0);
    P.dp=horzcat(P.dp,varargin{i}.dp);
    P.reg=vertcat(P.reg,varargin{i}.reg);
    P.Lx1=vertcat(P.Lx1,varargin{i}.Lx1);
    P.Lx2=vertcat(P.Lx2,varargin{i}.Lx2);
    P.dx=vertcat(P.dx,varargin{i}.dx); 
    P.a=horzcat(P.a,varargin{i}.a);
    P.b=horzcat(P.b,varargin{i}.b);
    
end
P.n=1:length(P.n);


end

