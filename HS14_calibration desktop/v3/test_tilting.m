% test the tilting
% here we found out that this actually works^^
im=tiltblock(:,:,1);


l1=600
l2=2000
dl=l2-l1+1
figure(50)
clf
subplot(1,10,1)

for i = 1: 10
    subplot(1,10,i)
    imbo4(im)
    
    ang=atand(f3_subpixcorr(im(:,l1),im(:,l2))/(dl));
    xlabel(num2str(ang))
    
    
    im=imrotate(im,-ang,'bilinear','crop');
end
    
