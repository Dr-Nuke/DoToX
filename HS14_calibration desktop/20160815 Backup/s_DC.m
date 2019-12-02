%% read ins the DC images.
% note: its the scatter - corrected DC images (big block of PE in front of
% the camera)

dc_empty_path      = strcat(readpath,'CD\');
dc_file='Boratedpoly_block_';
DC = zeros(length(xpix),length(ypix));
imax=5;
for i =1:imax,
    dc_filename=strcat(dc_empty_path,strcat(dc_file,num2str(i),'.fits')); 
    im=im2double(fitsread(dc_filename)');
    im=im(xpix,ypix);
    im=f_MediFilter1(im,thresh_medi);
    kernel=f_KernelGen(9,9,7.5);
    im= f_MediFilter2( im,kernel,0.06,'median' );    
    
    DC=DC+im
end
DC=DC/imax;