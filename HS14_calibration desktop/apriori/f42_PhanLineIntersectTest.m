%function [] = f42_PhanLineIntersectTest()

% test the PhanLineIntersect function
P=f42_ChannelPhanTransform;
%P=f4_ChannelPhan2;
P=f42_RotPhan(P,deg2rad(45));
fig=figure(67);
clf
set(fig,'Position',[43,500,840,420]);
subplot(1,3,1:2)

% subplot(1,2,1)
ax=gca;
% plot phantom
f42_PlotPhan(P,ax,[0,0,0],[1]);
% make & plot a line
sstart=[-1000,0];

n=101;
detpix=linspace(-15,15,n);

    c=[[0.7,1,0.7];
        [0.4,0.4,0.4];
        [0,0,1];
        [1,0,0.5]];
    
% myvalues forgammacross sections
% region swich: 1:air 2:AL 3:D2 4: inside (ChCl3vapor)
%in 1/mm
xs=[0,10.143,45.039,0.485]/100;


% the attenuation
MuS=zeros(1,n);
    
for i=1:n
    eend=[15,detpix(i)];
    l = f42_CreateLineObj(sstart,eend);



    %f42_PlotLineObj(l,ax,[1,2,3,4],[0,0,1])

    [t,xy,n_seg,regio,s] = f42_PhanLineIntersect(P,l);
    
    MuS(i)=exp(-xs(regio)*s);


    

    xy2=[sstart;xy;eend];



    for j=1:size(xy2)-1
        l=f42_CreateLineObj([xy2(j,1),xy2(j,2)], [xy2(j+1,1),xy2(j+1,2)]);
        f42_PlotLineObj(l,ax,[1],c(regio(j),:));
    end
end
    % end
%%
xlim([-15,15])
ylim([-15,15])

subplot(1,3,3)
ax2=gca()
barh(MuS)
ylim([1,n])
set(gca,'Position',[0.6 0.1100 0.2134 0.8150])
ax2.YTick=[]
