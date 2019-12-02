function [l] = f42_CreateLineObj(xy1,xy2)
% creates a line object from two points


l.xy1=[xy1(1),xy1(2)];
l.xy2=[xy2(1),xy2(2)];
l.dxy=xy2-xy1;
    
if (xy2(1)-xy1(1))~=0 %catch vertical line


    l.a= (xy2(2)-xy1(2))/(xy2(1)-xy1(1));
    l.b= xy1(1,2)-l.a*xy1(1,1);

else
    l.a=inf*sign(xy2(2)-xy1(2));
    l.b=0;
end

l.l=norm(l.xy2-l.xy1);

end

