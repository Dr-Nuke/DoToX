function [x,y] = f42_CheckAngle(x,y,c,p0,dp)
% checks if a given (intersection point lies within the angular section of
% a segment
% x = [x1,...,xn] of the questionable ponints. 
% y = [y1,...,yn] of the questionable ponints.
% c = [cx,cy] center point
% p0 = start angle
% dp = angular section size

ind=ones(length(x));

for i=1:length(x)
    a=atan2(y(i)-c(2),x(i)-c(1));
    p1=p0+dp;
    
    for ang=2*pi*[0,1,-1]   % for loop required to fix negative and
                            % >2pi angles
        % check if the point is in the right angular range
        if (min(p0,p1)<=a+ang) && (max(p0,p1)>=a+ang)
            ind(i)=0;
            continue
        end
    end
end

x=x(ind==0);
y=y(ind==0);

            
            
            
        


end

