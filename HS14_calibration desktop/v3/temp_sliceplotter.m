k=600 %slice number 600


% old centering
cent_old=shift(k)

smin=-5
smax=5
sn=11
extrashift=linspace(smin,smax,sn);
shift_test=cent_old+extrashift;


sino_orig=squeeze(block(:,k,:));

for i=1:sn
    sino=f3_center(sino_orig,shift_test(i));
    testrecon(:,:,i)=iradon(sino(:,1:end-1),angles,'spline','Hann',1,999);
    
    figure(i)
    imshow(squeeze(testrecon(:,:,i)),[])
    title(shift_test(i))
end
    
%%
figure()
for i=500:600
    imshow(squeeze(recontot(:,:,i)),[])
end

%%

i=1000;
sino=squeeze(block(:,i,:));
    sino=f3_center(sino,shift(i));
    sino=circshift(sino(:,1:end-1),11,2);
    

        %hardcode cropping
        sino=sino(201:end-200,:);
    
    testrecon=iradon(sino,angles,'spline','Hann',1,999);
    figure(21)
    imshow(testrecon,[])