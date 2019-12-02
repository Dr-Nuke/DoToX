function M = f4_DoseArea(M)
% does the beam non-uniformity correction

mask=zeros(size(M.im(:,:,1)),'single');
for i=1:size(M.dose,2)
    mask(M.dose(1,i):M.dose(2,i),M.dose(3,i):M.dose(4,i))=1;
end
M.mask = mask;

if M.proofs==1 % check the areas
    h=figure(4);clf;
    set(gcf,'name','dose area mask check')
    subplot(1,3,1)
    imshow(squeeze(M.im(:,:,1))',[]);set(gca,'YDir','normal')
    title(sprintf('%02d original',M.i))
    subplot(1,3,2)
    imshow(mask',[]);set(gca,'YDir','normal')
    title(sprintf('%02d mask',M.i))
    subplot(1,3,3)
    imshow((single(squeeze(M.im(:,:,1))).*(0.5+mask/2))',[]);set(gca,'YDir','normal')
    title(sprintf('%02d combined',M.i))
    fname=sprintf('Fig04 %02d dose mask',M.i);
    f_BoFig2PDF(h,fname)  
    
    %check if nothing rotates in these areas
    h=figure(5);clf;
    for i=1:M.nFrame
       a=single(squeeze(M.im(:,:,i))).*mask;
       dosemin(i)=min(a(a>0));
       dosemean(i)=mean(a(a>0));
    end
    
    subplot(2,1,1)
    plot(dosemean)
    title(sprintf('%02d mean dose in dose areas per frame',M.i));
    ylabel('mean dose')
    xlabel('frame number')
    grid on
    
    subplot(2,1,2)
    plot(dosemin)
    ylabel('minimum dose image value')
    xlabel('frame number')
    title(sprintf('if this has major dips, then some object \n rotates into the dose mask'))
    grid on
    fname=sprintf('Fig05 %02d dose area checks',M.i);
    f_BoFig2PDF(h,fname) 
    
    
end


end

