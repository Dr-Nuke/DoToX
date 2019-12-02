%Isosurface test sript

clear all
close all
clc
load('C:\data\Tomo_HS14\processed\recon.mat')
a=(size(recon));
x=a(1);
y=a(2);
z=a(3);
[X,Y,Z]=meshgrid(single(1:x),single(1:y),single(1:z));

ma=max(recon(:));
mi=min(recon(:));
%%
figure(1)
isovalue = 0.4*(ma-mi)+mi;
surf1 = isosurface(X,Y,Z,recon,isovalue);
p1 = patch(surf1);
isonormalsBo(x,y,z,recon,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(3); 
camlight; lighting gouraud

% isovalue = 0.4*(ma-mi)+mi;
% surf2=isosurface(x,y,z,recon,isovalue);
% p2 = patch(surf2);
% isonormals(x,y,z,recon,p2);
% set(p2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2);
% 
% isovalue = 0.6*(ma-mi)+mi;
% surf3=isosurface(x,y,z,recon,isovalue);
% p3 = patch(surf3);
% isonormals(x,y,z,recon,p3);
% set(p3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3);
% 
% isovalue = 0.8*(ma-mi)+mi;
% surf4=isosurface(x,y,z,recon,isovalue);
% p4 = patch(surf4);
% isonormals(x,y,z,recon,p4);
% set(p4,'FaceColor','blue','EdgeColor','none','FaceAlpha',1);
%----------
