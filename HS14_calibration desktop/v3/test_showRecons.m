fig=figure(18);


x=[350,2115];
y=[0.05,1]
xq=[500,1000,1500,2000];



writerObj = VideoWriter('YourAVI.avi');
clear M

for i=[1:10:900]
%imshow(squeeze(recontot(:,:,i)),[0,0.01]);

titlestring=sprintf('pixel plane %4d',i);

if and(350<i,i<2115)
    mm=interp1(x,y,i);
    titlestring=strcat(titlestring,sprintf(', LFT = %1.3fmm',mm));
end
imshow(squeeze(recontot(:,i,:)),[0,0.01]);
title(titlestring)
drawnow;

M(i)=getframe(fig);

end
%%

video=VideoWriter('flythrough','Motion JPEG 2000');
video.FrameRate=30;
open(video)
writeVideo(video,M);
close(video)
%%

ind_i=[1:50:2400]
for i=1:2
        sino=squeeze(block(:,ind_i(i),:)); %sinogram
        sino=f3_center(sino,shift(ind_i(i)));

        sinotest{i}=sino(201:end-200,:);
        try
            recontest(:,:,i)=iradon(sinotest{i}(:,1:end-1),angles,'spline','Hann',1,999);
        catch
            disp(sprintf('%4d went wrong',k))
        end
        %sinotot(:,:,k)=sino;
        end