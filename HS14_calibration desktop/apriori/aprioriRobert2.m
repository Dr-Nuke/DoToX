%make phantom

lft=0.2;% mm lft
res=0.01; % sampling of the bitmap pixels in mm

domain=30; %mm
nx=domain/res+1;  % number of x-pixels
ny=domain/res+1;  % number of y-pixels
% generate phantom, focus on 1 quarter
xcoord=linspace(-domain/2,domain/2,nx);
ycoord=linspace(-domain/2,domain/2,ny);
[xx,yy]=ndgrid(xcoord,ycoord); % coordinates matrix
P=nan(nx,ny); % phantom image

Ph = f42_ChannelPhanTransform_w_film(lft)
figure();clf;
f42_PlotPhan(Ph,gca,[1 0 0],[1]);

% 0 void
% 1 alu
% 2 h2o
% 3 chcl3 liq
% 4 chcle vap

P(1,1)=0;
for i=1:nx
    for j=1:ny
        
        
    end
end