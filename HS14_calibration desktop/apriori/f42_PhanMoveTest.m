function  []= f42_PhanMoveTest
% test the Phan Move function

P=f42_ChannelPhanTransform;

figure(1)
clf
ax1=gca;
f42_PlotPhan(P,ax1,[0,0,0],[1]);
hold on

P2=f42_PhanMove(P,5,0);
f42_PlotPhan(P2,ax1,[1,0,0],[1]);

P3=f42_PhanMove(P,0,5);
f42_PlotPhan(P3,ax1,[0,1,0],[1]);

P4=f42_PhanMove(P,-5,-5);
f42_PlotPhan(P4,ax1,[0,0,1],[1]);


end

