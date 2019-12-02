figure(30)
clf
gca()
hold on

figure(31)
clf
gca()
hold on

I=400;

imshow(squeeze(block(:,:,I))',[]);set(gca,'YDir','normal')
hold on
for i = I
    for j=3
        for k=1:3:F.n_angles 
            figure(31)
            plot([F.cfit(i,j,1),F.ProfEndPnt(i,j,k,1)],...
                 [F.cfit(i,j,2),F.ProfEndPnt(i,j,k,2)])
            figure(30)
            
            plot3(1:F.r_path+1,... % x
                k*ones(1,F.r_path+1),...% y
                squeeze(F.ImProfiles(i,j,k,:))) %z
        end
    end
end
    
figure(32)
clf
plot(squeeze(F.ImProfiles(i,j,1,:)))

%%


figure(35)
clf
gca();
hold on

series=80:420;
smax=length(series);

col=hsv(smax);
l=size(F.ImProfiles,4)
x=(1:l)/res;

for i=1:smax
    
    plot(x,squeeze(F.ImProfiles(series(i),1,70,:)),...
        'DisplayName',sprintf('%4.0f um',f_plane2lft(refac*series(i))),...
        'Color',col(i,:))
    
end
grid on
plot([5.14,5.14],[0,0.03],'Color',[0.5,0.5,0.5],'DisplayName','rod surface',...
    'LineWidth',2)

legend('Location','northwest')
xlabel('path distance, [mm]')
ylabel('image intensity')
title('sample curves for LFT calculation')

%%


figure(36)
clf
gca();
hold on


for i=1:smax
    y(i)=f_plane2lft(refac*series(i));
    xx=F.ImProfiles(series(i),:,:,:);
    LFTsum(i)=sum(xx(:));
        
    
end
scatter(LFTsum,y)%... %x value
    
grid on




    