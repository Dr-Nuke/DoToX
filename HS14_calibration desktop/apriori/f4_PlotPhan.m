function [] = f4_PlotPhan(phan,c)
% makes a plot of the phantom with correct segment heights
% phan = the phantom 
% c = color specifier, 
% fid = figure id


cmax=max(phan.c);
cmin=min(phan.c);
cmap=jet
fa=0.1; %face alpha


hold on
for i = 1:length(phan.c)
    col=f4_Colorgen(cmap,cmin,cmax,phan.c(i));
    
    %plot values
    plot3([phan.xy(1,i),phan.xy(1,i+1)],...
          [phan.xy(2,i),phan.xy(2,i+1)],...
          [phan.c(i),phan.c(i)],'Color',col,...
          'LineWidth',2);
      % add "bar faces"
      x=[phan.xy(1,i),phan.xy(1,i),phan.xy(1,i+1),phan.xy(1,i+1)];
      y=[phan.xy(2,i),phan.xy(2,i),phan.xy(2,i+1),phan.xy(2,i+1)];
      z=[0,phan.c(i),phan.c(i),0];
      
      fill3(x,y,z,'k','FaceAlpha',fa,'EdgeColor','none');
      
      
      ii=mod(i-2,length(phan.c))+1;
      %disp([i,ii])
      ql(i)=max(phan.c(i),phan.c(ii));

end
plot3(phan.xy(1,:),phan.xy(2,:),zeros(length(phan.xy(1,:))),'k')



grid on
xlabel('x')
ylabel('y')
zlabel('local attenuation')



end
