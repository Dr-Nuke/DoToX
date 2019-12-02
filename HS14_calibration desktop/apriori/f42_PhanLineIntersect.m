function [t,xy,n_seg,regio,s] = f42_PhanLineIntersect(P,l)
% this function returns the data on intersections between a Phatom and a
% given line
% P = Phantom
% l = line object

% t = vector containing the fractionals of the line distance where
% intersections appear
% xy = list of points of intersection
% n_seg = list of Phantom segment interceptions
% regio = list of material indicators
% s     = lengths of regions 

%initialize
t=[];
xy=[];
mat=[];

% 1) distinguish between line and arc

    for i=1:length(P.n)
        %disp(i);
        %disp(i)
        tt=[];
        xyt=[];
        matt=[];
        if P.CL(i)==1 %if the secment is an arc
            % arc intersection check
            [x,y]=linecirc(l.a,l.b,P.c(i,1),P.c(i,2),P.r(i));
            % [x,y] may describe now a list of 0-2 points
            
            % check if  and both points are identical (tangent case of
            % lincirc) or if no intersection point exitst
            if any(isnan(x)) % nan means linecirc returned no intersection
                x=[];
                y=[];
            elseif (x(1)== x(2)) && (y(1) == y(2)) % the tangent case
                % only keep the first point, second is identical
                x(2)=[];
                y(2)=[];
                
            else  % it found 2 intersections
               % do nothing
            end
            % remove points outside the angular section
            [x,y] = f42_CheckAngle(x,y,P.c(i,:),P.p0(i),P.dp(i));

            % find the according t, xy and mat (t for temporary)

            for j=1:length(x)
                
                ttn=(x(j)-l.xy1(1))/(l.xy2(1)-l.xy1(1));
                if (ttn>=0) && (ttn<=1)
                
                    tt=[tt,(x(j)-l.xy1(1))/(l.xy2(1)-l.xy1(1))];
                    xyt=[xyt;[x(j),y(j)]]; 
                    matt=[matt;P.reg(i,:)];
                end
                
            end
        
        % end of line-cicrle interaction
        else % else if it is a straight line
            % line intersection check
            [tt,~]=f4_intercept(l.xy1,l.xy2,P.Lx1(i,:),P.Lx2(i,:));
            
            %check if there is an intersection <=> 0<=tt<=1
            if all([tt(1)>=0,tt(1)<=1,tt(2)>=0,tt(2)<=1])
                tt=tt(1);
                matt=[P.reg(i,:)];
                xyt=l.xy1+tt*l.dxy;
            else
               tt=[];
            end


        end
    
    % now append t,xy,nseg and mat 
    t=[t;tt(:)];
    xy=[xy;xyt];
    mat=[mat;matt];


% at last sort the outputs for increasing t

% re-produce the materials
   
    
    end
    if length(t)>0
        [t,tsort]=sort(t); %sort intersections in beamdirection
        %xy=xy(tsort);
        xy=[xy(tsort,1),xy(tsort,2)];
        mat=[mat(tsort,1),mat(tsort,2)];
    end
    n_seg=length(t);
    try
    regio=f42_RegioFinder(mat,1); % the regionswiches
    catch
        disp(mat)
    end
    s=([t;1]-[0;t])*l.l;     % the length of each region
    
end

