% compare a measurement ile with a calculated file

load('d2osin')

a=MuS;
b=d2osin;
f1=figure(1);
set(f1,'Name','raw data')
clf
subplot(2,1,1)
imshow(a,[])
title('simulation')
subplot(2,1,2)
imshow(b,[])
title('hs14 data')

b=b(:,330:1030);

f2=figure(2);
set(f2,'Name','cropped data')
clf
subplot(2,1,1)
imshow(a,[])
title('simulation')
subplot(2,1,2)
imshow(b,[])
title('hs14 data')

%% norm
f3=figure(3);
set(f3,'Name','normalized data')
clf

a=(a-min(a(:)))/(max(a(:))-min(a(:)));
b=(b-min(b(:)))/(max(b(:))-min(b(:)));
subplot(2,1,1)
imshow(a,[])
title('simulation')
subplot(2,1,2)
imshow(b,[])
title('hs14 data')


%% Histogram
f4=figure(4);
set(f4,'Name','histograms')
clf

nhist=100;
subplot(2,1,1)
hist(a(:),nhist);
title('simulation')
subplot(2,1,2)
hist(b(:),nhist);
title('hs14 data')
%%

B = imresize(MuS,size(d2osin));

f5=figure()


















