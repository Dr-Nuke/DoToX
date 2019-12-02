xrang=1:751;


zrang=round(linspace(1,1225,751));

figure(1)
clf
for i=1:751
    disp(i)
    
    
    %subplot(1,3,[1,2])
    imshow(squeeze(recontot2(:,:,zrang(i))),[]);
    title(i)
    drawnow
%     subplot(1,3,3)
%     imshow(squeeze(recontot2(i,:,:)),[]);
    
%     M(i)=getframe(fig);
   %     video=VideoWriter(['C:\Users\robersl\Downloads\video_3.mp4'],'MPEG-4');
%     video.FrameRate=10;
%     open(video)
%     writeVideo(video,M);
%     close(video)

end