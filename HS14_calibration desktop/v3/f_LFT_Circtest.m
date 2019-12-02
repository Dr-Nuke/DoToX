function [ output_args ] = f_LFT_Circtest(block,k,F)
im=squeeze(block(:,:,k));
%im(1:10,1:30)=max(im(:)); %orientation test   
imbw=im2bw(f_normalize(im),0.6);
[c,rad,m]=f_hugh_4(imbw,F.r_min,F.dr,F.S,F.o,0.6);


figure(17)
clf

imshow(im',[]) ;
set(gca,'YDir','normal')
title('LFT circtest')
hold on

%add centers
plot(c(:,1),c(:,2),'xr')

%add circles



end

