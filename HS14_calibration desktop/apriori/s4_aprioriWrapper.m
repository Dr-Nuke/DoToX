
close all
clear all
clc
phan=f4_ChannelPhan(111);
iimin=1;
iimax=1;
for ii=[iimin:iimax]%1:11
    disp(ii)
    figid=ii;
    noiselvl=0.3;  % 0... 0.5
    geom.gamma=10*ii;
    apriori
end

%%

figure(55)
clf
subplot(3,3,[1:6])

plot3(phan.xy(1,1:end-1),phan.xy(2,1:end-1),phan.c,'DisplayName','phantom')

grid on
xlabel('x')
ylabel('y')
zlabel('local attenuation')
legend('show')
title(sprintf('%d phantom segments, values randomly in [1,2]',length(phan.c) ))

%title(sprintf('%d pixels, %d angles, %d object segments',...
%    geom.nd,geom.gamma,length(phan.c)))

%%
figure(56)
clf

c=jet(11);
for ii=iimin:iimax
% log plot

%semilogy(abs(((xx(i,:)-phan.c)./phan.c)'),'DisplayName',num2str((i-1)/20),'Color',c(i,:))
semilogy(abs(((xx(ii,:)-phan.c)./phan.c)'),'DisplayName',num2str(10*ii),'Color',c(ii,:))
hold on
drawnow()
end
grid on

xlabel('object segment#')
ylabel('relative deviation ')
xlim([1,max(length(x),2)])
title({'abs(relative error) for different detector pixel numbers,';...
       'detector size 40mm'})

% title({'abs(relative error) for different noise levels,';...
%        'with phantom values randomly in [1,2]'})
legend('show')


%%
figure(46)
clf

c=jet(11);
for ii=iimin:iimax
% log plot

semilogy(abs(((xx(ii,:)-phan.c)./phan.c)'),'DisplayName',num2str((ii-1)/20),'Color',c(ii,:))
hold on
drawnow()
end
grid on

xlabel('object segment#')
ylabel('relative deviation ')
xlim([1,max(length(x),2)])
title({'abs(relative error) for different detector pixel numbers,';...
       'detector size 40mm'})

% title({'abs(relative error) for different noise levels,';...
%        'with phantom values randomly in [1,2]'})
legend('show')

