function [ output_args ] = f4_PlotResult(phan,x,fid,geom,noise)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
figure(fid)
clf

x=reshape(x,length(x),1);
phan.c=reshape(phan.c,length(x),1);
dx=(x-phan.c)./phan.c;

%% 3d plot
subplot(4,3,[1:6])

% plot3(phan.xy(1,1:end-1),phan.xy(2,1:end-1),phan.c,'DisplayName','phantom') %old
% plot3(phan.xy(1,1:end-1),phan.xy(2,1:end-1),x,'x','DisplayName','result') % old 

cmax=max(abs(dx));
cmin=-cmax;
cmap=[flipud(jet);jet];%colormap([flipud(jet);jet]);
fa=0.1; %face alpha
hold on
for i = 1:length(phan.c)
    %plot phantom values
    p1=plot3([phan.xy(1,i),phan.xy(1,i+1)],...
          [phan.xy(2,i),phan.xy(2,i+1)],...
          [phan.c(i),phan.c(i)],'k:',...
          'LineWidth',2);
      % add "bar faces"
      xx=[phan.xy(1,i),phan.xy(1,i),phan.xy(1,i+1),phan.xy(1,i+1)];
      yy=[phan.xy(2,i),phan.xy(2,i),phan.xy(2,i+1),phan.xy(2,i+1)];
      zz=[0,phan.c(i),phan.c(i),0];
      
      fill3(xx,yy,zz,'k','FaceAlpha',fa,'EdgeColor','none');
      % add the calculated values
      p2=plot3([phan.xy(1,i),phan.xy(1,i+1)],...
          [phan.xy(2,i),phan.xy(2,i+1)],...
          [x(i),x(i)],'Color',f4_Colorgen(cmap,cmin,cmax,dx(i)),...
          'LineWidth',2);      

end
colormap(cmap);
c=colorbar();
c.Label.String = 'the color of the a priori segments indicate their relative error';
caxis([cmin,cmax]);
% add the ground contour
plot3(phan.xy(1,:),phan.xy(2,:),zeros(length(phan.xy(1,:))),'k')

grid on
xlabel('x')
ylabel('y')
zlabel('local attenuation')
legend([p1,p2],'phantom','a priori computation')
title(sprintf('%d pixels, %d angles, %d object segments, added noise lvl %.2f',...
    geom.nd,geom.gamma,length(phan.c),noise))


%% log plot
subplot(4,3,7:9)
semilogy(abs((x-phan.c)./phan.c),'DisplayName','|deviation|')
hold on
grid on

xlabel('object segment#')
ylabel('relative deviation ')
xlim([1,max(length(x),2)])

legend('show')

%% normal plot
subplot(4,3,10:12)
plot((x-phan.c)./phan.c,'DisplayName','deviation')
hold on
grid on

xlabel('object segment#')
ylabel('absolute deviation')
xlim([1,max(length(x),2)])
legend('show')


end

