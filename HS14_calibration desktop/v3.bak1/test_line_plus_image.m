%% plot line profiles
figure(2)
clf
imagesc(cl33');
set(gca,'YDir','normal'); 
colormap('gray')
xlabel('x');
ylabel('y');
zlabel('z');
col=hsv(n);
 

cd = [uint8(jet(n)*255) uint8(ones(n,1))].' ;%'
for i=1:4



    hold on
    for j=1:n
        scatter3(squeeze(cx_tot(i,j,:)),...
              squeeze(cy_tot(i,j,:)),...
              squeeze(c_tot(i,j,:)),...
              2,...
              squeeze(c_tot(i,j,:)),'.');
          colormap(hsv)
        
    end
        


end




