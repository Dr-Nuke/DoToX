


for k=100%:100:900

im=squeeze(block(:,:,k));  

imbw=im2bw(f_normalize(im),0.4);
[c(k,:,:),rad(k,:),m(k,:)]=f_hugh_4_2(imbw,40,5,F.S,F.o,0.4); %have this file in the same folder
            

    
    figure(35)
    clf
    imshow(squeeze(im'),[]);set(gca,'YDir','normal') % image
    hold on
    %scatter(F.c(k,:,1),F.c(k,:,2),'gx') % center points
    scatter(F.c(k,:,1),F.c(k,:,2),'rx') 
    plot(c(k,:,1),c(k,:,2),'b')
    plot(F.cfit(k,:,1),F.cfit(k,:,2),'g')
    title(k)
    %waitforbuttonpress;
    pause(0.1)
    figure(32);imshow(imbw,[])

end

%%
    figure(35)
    viscircles(c,rad)

    
    %%
[c,rad,m] = imfindcircles(im2,[40,45],'Sensitivity',0.99,'ObjectPolarity','bright')
figure(36);clf;
imshow(im2',[]);set(gca,'YDir','normal')
viscircles(fliplr(c),rad)