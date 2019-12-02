fh=figure(312);clf
imshow(mtomo(:,:,1),[])
set(fh,'units','normalized',...
    'outerposition',[0 0 1 1])

hold on

colorbar
axis on
P=f42_ChannelPhanTransform();
P2=f42_RotPhan(P,deg2rad(43.7));
P2=f42_PhanZoom(P2,7.9);
P2=f42_PhanMove(P2,134.5,134.4);
f42_PlotPhan(P2,gca,[1 0 0],[1 4])

%%
for i =[62 46 29 13];1:length(P2.c)
    plot(P2.c(i,1),P2.c(i,2),'ob');
    text(P2.c(i,1),P2.c(i,2),num2str(i),'FontSize',20)
end
