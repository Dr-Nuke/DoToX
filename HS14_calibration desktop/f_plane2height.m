function height= f_plane2height(plane,cut)
% gives the HS14 plug heigt in cm based on pixel plane    
% plane = value or array of pixel plane values
% cut = cut flag. 1: cut off at top & bottm. 0: just use linear model


s=size(plane); %could be a vector

x=[354 2127]; % manual guess of the start and end pixel planes of the plug

y=[0 8]; % height of the plug
a=(y(2)-y(1))/(x(2)-x(1)); % inclination


b=y(1)*ones(s)-a*x(1)*ones(s);
height=a*plane+b;

if cut
    height(height>8)=8;
    height(height<0)=0;
end
% if plane<x(1)
%     height=0
% elseif plane > x(2)
%     height=8
% else


    
    
    



end



