load('G:\cbolesch\for robert\equifilm.mat');
x=equifilm;
a = min(x(:));
b = max(x(:));

figure()
for i=1:3
ax(i)=subplot(1,3,i);imshow(x(:,:,i),[a b]);
    
    
end

linkaxes(ax)