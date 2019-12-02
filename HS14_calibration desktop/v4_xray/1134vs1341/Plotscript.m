%% plotscriiiiiipt!

%load('LFT HS14 diff.mat')

%% 1 plot of the LFT map

figure(40);clf;
text(0.5,0.5,'test')
xtx=round(linspace(1,F.n_angles,7));
xlab={'','-45','',0,'','45',''};

for i=1:4
subplot(1,4,i)
imshow(squeeze(F.LFT(:,i,:)),[2 10]/1000);
set(gca,'YDir','normal')
title(sprintf('rod %u',i))
xlabel(sprintf('ang. pos.\n [deg]'))
    axis on
    set(gca,'ytick',[])
    
xticks(xtx)
xticklabels(xlab)    
end

% add one y axis labelz
subplot(1,4,1)
ylabel('rod height, [pixel]')
%ytx=round(f_PixReCalc(linspace(354,2127,5),refac));
for i=1:5
    ylab{i}=sprintf('%.0f',f_plane2height(ytx(i)*refac,1));

    
end
%yticks(ytx)
%yticklabels(ylab)

% overall title



%% LFT axial profile plot per rod
figure(41);clf;gca();
hold on

x=1:imrange(2,2);
xlft=f_plane2lft(x*refac);
xheight=f_plane2height(x*refac,0);


for i=1:4
    plotlft=mean(F.LFT,3);
    plot(plotlft(:,i),1:size(F.LFT,1),'displayname',...
        sprintf('Rod %d',i))
end
grid on
xlbl='mean attenuation integral';
xlabel(xlbl)
title(sprintf('mean LFT signal per rod vs rod height \n X-Ray annular flow'))
ylabel('rod height [pixel]')
legend()
%%
figure('Name', 'median filter for 4 rods');clf;gca();hold on;

yy=mean(F.LFT,3);
yyy=medfilt1(yy);
y=yy;
filt=(abs(yyy-y)./y)>0.05;
y(filt)=yyy(filt);
disp([i,sum(filt)])



for i=1:4
    subplot(2,2,i)
    hold on
    plot(yy(:,i))
    plot(y(:,i))
    title(sprintf('rod %u',i))
    xlim([0 500])
    %ylim([-0.1,2])
    
    if or(i==2,i==4)
        set(gca,'yticklabel',[])
        ylabel(ylbl)
    end
    if or(i==1,i==2)
        set(gca,'xticklabel',[])
        xlabel('plug height')
    end
    grid on 
    
end
grid on

cl3crop=180:420; %manually extracted from the graph
d2ocrop=97:156;


%% LFT axial profile plot per rod, with median filter
figure(41);clf;gca();
hold on

x=1:imrange(2,2);
xheight=f_plane2height(x*refac,0);
xlft=f_plane2lft(x*refac);


for i=1:4
    plot(xheight,y(:,i));
end
grid on
ylbl='mean attenuation integral';
ylabel(ylbl)
title('blubb')

%xlim(round([250,2200]/refac))
%ylim([-0.5,1.5])


xfitw=[450,775];
xfitwrange=f_plane2height(xfitw(1):xfitw(2),0); % WATER FIT

xfitc=[950,2100];
xfitcrange=f_plane2height(xfitc(1):xfitc(2),0); % WATER FIT

line([xfitwrange(1),xfitwrange(1)],[0,0.2],'color','b')
line([xfitwrange(end),xfitwrange(end)],[0,0.4],'color','b')

line([xfitcrange(1),xfitcrange(1)],[0,0.5],'color','r') 
line([xfitcrange(end),xfitcrange(end)],[0,1.1],'color','r')

xfitwRange2=round(xfitw(1)/refac):round(xfitw(2)/refac);
xfitcRange2=round(xfitc(1)/refac):round(xfitc(2)/refac);

xnew=linspace(0,8,101);

for i =1:4
    pw(i,:)=polyfit(xheight(xfitwRange2),y(xfitwRange2,i)',1);
    pc(i,:)=polyfit(xheight(xfitcRange2),y(xfitcRange2,i)',1);
    
    yfitw(i,:)=polyval(pw(i,:),xnew);
    yfitc(i,:)=polyval(pc(i,:),xnew);
        
    %plot(xnew,yfitw(i,:))
    %plot(xnew,yfitc(i,:))
    
end

    
ylim([-0.2 1.2])
xlabel('Rod height [cm]')









