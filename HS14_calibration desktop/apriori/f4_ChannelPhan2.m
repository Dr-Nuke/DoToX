function [P] = f4_ChannelPhan2
% creates the second generation Channel phantom
    [p1,p3,x2,y2] = f42_GeomSolver(); % obtain the numbers for the second point

    %% geometric data from inventor CAD
                           % radii
    r1=13.604;
    r2=2;
    r3=5.14;
    
    
    p2=p1;
    
    dp1=p1;
    dp2=p3-pi-p1;
    dp3=-(p3-(5/4)*pi);
    
    
    %% create the geometry, first 1/8th, inner boundary

    P.n=[1];            % segment number
    P.CL=[1];            % Circle- Line indicator: equals 1 if the segment is a circle   
    P.c=[x2,y2];       % circle cter
    P.r=[r2];           % radius
    P.p0=[p2];          % circle segment start angle
    P.dp=[dp2] ;        % circle segment angular size
    P.reg=[4,2];        % region swich: 1:air (outside) 2:AL 3:D2 4: inside

    P.Lx1=[nan,nan]; %starting points
    P.Lx2=[nan,nan]; % ending points
    P.dx=[nan,nan]; % differential vector
    P.a=[nan];          % slope
    P.b=[nan];          % y-intercept
    
    
    
    % P2=hold on
    
    %% d2o channel boundary
    
    [p4 x4 y4 t] = f42_GeomSolver2(); % get the data of the 4th point
    
    x6=6.7+1/sqrt(2);
    y6=6.7+1/sqrt(2);
    dp4=pi*(2+1/4)-p4;
    x42=x6+t/sqrt(2);
    y42=y6-t/sqrt(2);
    dp32=p4-5/4*pi;
    
%     Q.n=[1,2,3];            % segment number
%     Q.CL=[1,1,0];            % Circle- Line indicator: equals 1 if the segment is a circle   
%     Q.c=[6.7,6.7;x4,y4;0,0];       % circle cter
%     Q.r=[4.14,2,0];           % radius
%     Q.p0=[5/4*pi,p4,0];          % circle segment start angle
%     Q.dp=[dp32,dp4,0] ;        % circle segment angular size
%     Q.reg=[3,2;...                     % region swich: 1:air 2:AL 3:D2O 4:inside  
%          3,2;...
%          3,2;...
%          ];
% 
%     Q.Lx1=[0,0;0,0;x42,y42]; %starting points
%     Q.Lx2=[0,0;0,0;x6,y6]; % ending points

    Q.n=[2];            % segment number
    Q.CL=[1];            % Circle- Line indicator: equals 1 if the segment is a circle   
    Q.c=[x4,y4];       % circle cter
    Q.r=[2];           % radius
    Q.p0=[p4];          % circle segment start angle
    Q.dp=[dp4] ;        % circle segment angular size
    Q.reg=[3,2];                    % region swich: 1:air 2:AL 3:D2O 4:inside  
        

    Q.Lx1=[nan,nan]; %starting points
    Q.Lx2=[nan,nan]; % ending points
    Q.dx=[nan,nan;nan,nan]; % differential vector
    Q.a=[nan,nan];          % slope
    Q.b=[nan,nan];          % y-intercept
    
P=f42_AddPhan(P,Q);   
    %% outter wall
    
    [x5 y5 p5 p45] = f42_GeomSolver3(x2,y2,x4,y4); % numbers for point 5
    

    
    dp25=p5-pi-dp1;
    dp5=-(p5-(p45+pi));
    dp45=pi/4-p45;


    R.n=[2,3,4];            % segment number
    R.CL=[1,1,1];            % Circle- Line indicator: equals 1 if the segment is a circle   
    R.c=[x2,y2;x5,y5;x4,y4];       % circle cter
    R.r=[3,2,3];           % radius
    R.p0=[p1,p5,p45];          % circle segment start angle
    R.dp=[dp25,dp5,dp45] ;        % circle segment angular size
    R.reg=[1,2;1,2;1,2];                    % region swich: 1:air 2:AL 3:D2O 4:inside  
         
    R.Lx1=[nan,nan;nan,nan;nan,nan]; %starting points
    R.Lx2=[nan,nan;nan,nan;nan,nan]; % ending points
    R.dx=[nan,nan;nan,nan]; % differential vector
    R.a=[nan,nan];          % slope
    R.b=[nan,nan];          % y-intercept    
    
    
P=f42_AddPhan(P,R);      

P2=f42_PhanMirr45(P);

% now adding the segments that are not suited for 1/8th unfolding 
% i.e. extend over n*45° lines

    % first the arcs

    dp3e=-2*(p3-5/4*pi);

    S.n=[1,2,3,4];            % segment number
    S.CL=[1,1,1,1];            % Circle- Line indicator: equals 1 if the segment is a circle   
    S.c=[0,0;0,0;6.7,6.7;6.7,6.7];       % circle cter
    S.r=[r1,r1+1,r3,r3-1];           % radius
    S.p0=[-p1,-p1,p3,p4];          % circle segment start angle
    S.dp=[2*p1,2*p1,dp3e,-2*dp32] ;        % circle segment angular size
    S.reg=[4,2;2,1;4,2;2,3];               % region swich: 1:air 2:AL 3:D2O 4:inside  
        

    S.Lx1=[nan,nan;nan,nan;nan,nan;nan,nan]; %starting points
    S.Lx2=[nan,nan;nan,nan;nan,nan;nan,nan]; % ending points
    S.dx=[nan,nan;nan,nan;nan,nan;nan,nan]; % differential vector
    S.a=[nan,nan,nan,nan];          % slope
    S.b=[nan,nan,nan,nan];          % y-intercept
    
    %then the line segments
    
        
    x10=6.7+1/sqrt(2)+t/sqrt(2);
    y10=6.7+1/sqrt(2)-t/sqrt(2);

    
    x11=6.7+1/sqrt(2)-t/sqrt(2);
    y11=6.7+1/sqrt(2)+t/sqrt(2);
    
    a10=(y11-y10)/(x11-x10);
    b10=y10-a10*x10;
    
    x12=6.7+2/sqrt(2)+t/sqrt(2);
    y12=6.7+2/sqrt(2)-t/sqrt(2);

    x13=6.7+2/sqrt(2)-t/sqrt(2);
    y13=6.7+2/sqrt(2)+t/sqrt(2);
   
    a12=(y13-y12)/(x13-x12);
    b12=y12-a12*x12;
    
    
    T.n=[1,1];            % segment number
    T.CL=[0,0];            % Circle- Line indicator: equals 1 if the segment is a circle   
    T.c=[nan,nan;nan,nan];       % circle cter
    T.r=[nan,nan];           % radius
    T.p0=[nan,nan];          % circle segment start angle
    T.dp=[nan,nan] ;        % circle segment angular size
    T.reg=[3,2;2,1];               % region swich: 1:air 2:AL 3:D2O 4:inside  
        

    T.Lx1=[x10,y10;x12,y12]; %starting points
    T.Lx2=[x11,y11;x13,y13]; % ending points
    T.dx=[x11-x10,y11-y10;x13-x12,y13-y12]; % differential vector
    T.a=[a10,a12];          % slope
    T.b=[b10,b12];          % y-intercept
    
    


P=f42_AddPhan(P,P2,S,T);
end



























