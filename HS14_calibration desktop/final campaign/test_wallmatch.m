% we try to match two reconstructions by the wall perimeter

% 1) get 2 absolute tomos
% 2) create a mask based on the gradient
% 3) optimize for alpha, x0,y0, zoom

% 1
%load('T')
%recon2=f.BoLoad(,T)


imref=squeeze(recon2(:,:,100,1));
im=squeeze(recon2(:,:,100,4));

imshow(im-imref,[]);

[Gmag,~] = imgradient(im);

imshow(Gmag,[]);