function lft= f_plane2lft(plane)
    

x=[354 2127];

y=[50 1000];

if plane<x(1)
    lft=0
elseif plane > x(2)
    lft=0
else


    a=(y(2)-y(1))/(x(2)-x(1));
    b=y(1)-a*x(1);

    lft=a*plane+b;
end

;;;
end



