function [] = f42_CheckAngleTest()
% testing the CheckAngle function
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['testing ',name,'...'])             %command line output  
    
    % range of random points

    
    % make random oints
    xmin=-1;
    xmax=1;
    ymin=-1;
    ymax=1;
    n=1000;
    x=(rand(1,n)-0.5)*(xmax-xmin);
    y=(rand(1,n)-0.5)*(ymax-ymin);
    
    % the empty case
    x=[];
    y=[];
    
    % center point
    c0=[1,1];

    figure(18);
    clf
    scatter(x,y,'rx');
    hold on
    
    nrot=101;
    ang=linspace(0,2*pi,nrot);
    dp=pi/4;
    
    for i=1:nrot
        disp(i)
        p0=ang(i);
        % forward angle
        [x2,y2] = f42_CheckAngle(x,y,c0,p0,dp);
        s=scatter(x2,y2,'b+');
    
        c1=c0+2*[cos(p0),sin(p0)];
        c2=c0+2*[cos(p0+dp),sin(p0+dp)];
        l1=line([c0(1),c1(1)],[c0(2),c1(2)]);
        l2=line([c0(1),c2(1)],[c0(2),c2(2)]);
        
        % and backward angle
        [x2,y2] = f42_CheckAngle(x,y,c0,p0,-dp);
        s2=scatter(x2,y2,'g+');
    
        c1=c0+2*[cos(p0),sin(p0)];
        c3=c0+2*[cos(p0+-dp),sin(p0+-dp)];
        %l1=line([c0(1),c1(1)],[c0(2),c1(2)]);
        l3=line([c0(1),c3(1)],[c0(2),c3(2)]);
        
        
        grid on
        axis equal    
        axis([-2,2,-2,2]);
        drawnow()
        pause(0.1);
        delete([s,s2,l1,l2,l3]);
        
    end
    
end

