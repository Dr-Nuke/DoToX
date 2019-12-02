function [ phan ] = f4_CircPhan(r,ds)
% makes a circle phantom.
% r = radius in mm
% n_circle = number of concentric circles
% ds= segment length





    u=2*pi*r;       % umfang/ circumference
    n_seg = floor(u/ds);   % anzahl segmente
    phi=linspace(0,2*pi,n_seg); %list of angles
    
    x=r*cos(1*phi);
    y=r*sin(1*phi);
    phan.xy=[x;y];

    phan.c=rand(1,n_seg-1)+1;
    
    
    
    %% hardcode
%     phan.c(1:4:end)=2;
%     phan.c(2:4:end)=2;

    
    %phan.c=cos(5*phi(1:end-1));
    

end

