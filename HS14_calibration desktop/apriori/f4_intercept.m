function [t,cosphi] = f4_intercept(xy10,xy11,xy20,xy21)
% finds the interception of two lines, defined by their end point coordinates 
% xy10 = [x,y] of the 1st point of the 1st line
% xy11 = [x,y] of the 2nd point of the 1st line
% xy20 = [x,y] of the 1st point of the 2nd line
% xy21 = [x,y] of the 2nd point of the 2nd line


vec1=[xy11(1)-xy10(1);... % vector of first line
      xy11(2)-xy10(2)];
  
vec2=[xy21(1)-xy20(1);... % vector of second line
      xy21(2)-xy20(2)];  
    

A=[vec1(1),-(vec2(1));... % lin.system matrix
   vec1(2),-(vec2(2))];

b=[-(xy10(1)-xy20(1)),... % lin.system result
   -(xy10(2)-xy20(2))];

t=A\b';                     % unknown t

% cos(phi)= a * b / (|a|*|b|)
cosphi=dot(vec1,vec2)/(norm(vec1)*norm(vec2));




end

