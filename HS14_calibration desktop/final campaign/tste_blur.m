img = im2double(blurred);
figure(1);clf;
subplot(1,2,1);
imshow(img), title('Original Image');



deblurSigma = 1; %Adjust this to get the most visually pleasing results

motion_noise = fspecial('gaussian', 15,deblurSigma);

luc1 = deconvlucy(img,motion_noise);

subplot(1,2,2); imshow(luc1);
title('Disk and Lucy');