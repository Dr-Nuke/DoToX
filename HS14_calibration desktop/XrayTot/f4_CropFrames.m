function M = f4_CropFrames(M,i)
% crops the data frames accordingly


if M.proofs
    h=figure(1);clf;

    a=squeeze(M.raw(:,45,:));
    imshow(a,[]);
    title(sprintf('%02d Raw sinogram of the angle gauge region',M.i))
    fname=sprintf('Fig01 %02d raw sinogram',i);
    f_BoFig2PDF(h,fname)

    h=figure(2);clf;
    imshow(squeeze(M.raw(:,45,M.d.crop1(i,1):M.d.crop1(i,2))),[])
    title(sprintf('%02d cropped block',i))
    fname=sprintf('Fig02 %02d raw cropped sinogram',i);
    f_BoFig2PDF(h,fname)    
    
end

M.im=M.raw(:,:,M.d.crop1(i,1):M.d.crop1(i,2));
M.nFrame(i)=size(M.im,3);
end

