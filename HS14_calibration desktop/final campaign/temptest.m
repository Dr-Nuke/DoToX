
figure(1347);clf;
cas=2;
pla=10;
pin=1;
%imshow(mtomo(:,:,cas,pla),[]);

[xx,yy]=ndgrid(1:269,1:269);
    zz=zeros(269,269);
su=surf(xx,yy,zz,mtomo(:,:,cas,pla),'edgecolor','none');
colormap(gray)
angs=1:10:F.n_angles;
nang=length(angs);
col=hsv(nang*4);
hold on
for pin=1:4
for ang2=1:nang;
    ang=angs(ang2);
    x=linspace(F.cen.cen(pin,1),F.Pa.ProfEndPnt(pin,ang,1),size(F.Profs,5));
    y=linspace(F.cen.cen(pin,2),F.Pa.ProfEndPnt(pin,ang,2),size(F.Profs,5));
    plot3(x,y,squeeze(F.Profs(cas,pla,pin,ang,:)),'color',col(ang2+(pin-1)*nang,:))

end
end