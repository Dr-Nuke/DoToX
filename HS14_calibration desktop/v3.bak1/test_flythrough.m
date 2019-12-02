% try to make a flythrough video

%given variable 'recon'

a=size(recon);

fig=figure(87);
clf
for i=1:a(3)
    imshow(squeeze(recon(:,:,i)),[0,0.01])
    disp(i)
    
    
    M(i)=getframe(fig);
end


video=VideoWriter(['C:\data\Tomo_HS14\flythrough.mp4'],'MPEG-4');
video.FrameRate=25;
open(video)
writeVideo(video,M);
close(video)
