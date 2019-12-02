ind_d=[1,2,3,4];

figure(1); clf; 
ax1=subplot(2,2,1);hold on;

P7=f42_ChannelPhanTransform();
f42_PlotPhan(P7,ax1,[1,0,1],ind_d);
%%

ax3=subplot(2,2,3); hold on;
P8=f4_ChannelPhan2();
f42_PlotPhan(P8,ax3,[0,0,0],ind_d);

%%
ax2=subplot(2,2,2); hold on;
c=hsv(4);
ang=linspace(0,2*pi,5);

%figure(2);clf; ax2=gca;

for i=[1,2,3,4]
    Prot=f42_RotPhan(P8,ang(i));
    f42_PlotPhan(Prot,ax2,c(i,:),ind_d);

end
    %%



% 
% P5=f42_RotPhan(P4,pi)
% f42_PlotPhan(P5,ax1,[0,1,1]);