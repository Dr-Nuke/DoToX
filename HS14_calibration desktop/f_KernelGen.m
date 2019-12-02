function [kernel] = f_KernelGen(d,do,di)
if mod(d,2)~=1;
    disp('bad kernel diameter not odd');
end
%pixel diameter of kernel . inner and outer. outer must be odd
ro=do/2; %outer radius of kernel
ri=di/2; %inner radius of kernel
r=((d+1)/2); % middle pixel
[x,y]=meshgrid(1:d,1:d);
kernel=zeros(d);
kernel(and((x-r).^2+(y-r).^2>=ri^2,(x-r).^2+(y-r).^2<=ro^2))=1;
end