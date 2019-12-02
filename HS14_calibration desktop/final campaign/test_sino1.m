clear sino1
load sino1

tstr={'orig','lr','ud','ud-lr','lr-ud','+90','-90'}
sino1(:,:,2)=fliplr(sino1(:,:,1));
sino1(:,:,3)=flipud(sino1(:,:,1));
sino1(:,:,4)=fliplr(flipud(sino1(:,:,1)));
sino1(:,:,5)=flipud(fliplr(sino1(:,:,1)));
sino2(:,:,1) = imrotate(sino1(:,:,1),90);
sino2(:,:,2) = imrotate(sino1(:,:,1),-90);




figure(1);clf;
imshow(a.FBPexplFan(sino1(:,:,1),319,ang,detPitch,src,det),[])
title('orig')

figure(2);clf;
imshow(a.FBPexplFan(sino1(:,:,2),319,ang,detPitch,src,det),[])
title('lr')

figure(3);clf;
imshow(a.FBPexplFan(sino1(:,:,3),319,ang,detPitch,src,det),[])
title('ud')

figure(4);clf;
imshow(a.FBPexplFan(sino1(:,:,4),319,ang,detPitch,src,det),[])
title('ud-lr')
%%
figure(5);clf;
imshow(a.FBPexplFan(sino1(:,:,5),319,ang,detPitch,src,det),[])
title('lr-ud')

figure(6);clf;
imshow(a.FBPexplFan(sino2(:,:,1)',319,ang,detPitch,src,det),[])
title('90')

figure(7);clf;
imshow(a.FBPexplFan(sino2(:,:,2)',319,ang,detPitch,src,det),[])
title('-90')


%%
c=2;
p=6;
r=1
detPitch=0.127;
det=50;
src=950;
sino=f.fraccircshift(squeeze(flip(tra(:,p,:,c,1),3)),...
            -shift_raw(p,c,r));
        figure(23);clf;
imshow(a.FBPexplFan(sino',recsize,ang,detPitch,src,det),[])


        
        
































