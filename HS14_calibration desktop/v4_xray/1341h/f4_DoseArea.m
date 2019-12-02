function M = f4_DoseArea(M)
% does the beam non-uniformity correction

mask=zeros(size(M.im(:,:,1)),'single');
for i=1:size(M.dose,2)
    mask(M.dose(1,i):M.dose(2,i),M.dose(3,i):M.dose(4,i))=1;
end
M.mask = mask;

if M.proofs==1 % check the areas
    figure(4);clf;
    set(gcf,'name','dose area mask check')
    subplot(1,3,1)
    imshow(squeeze(M.im(:,:,1))',[]);set(gca,'YDir','normal')
    title('original')
    subplot(1,3,2)
    imshow(mask',[]);set(gca,'YDir','normal')
    title('mask')
    subplot(1,3,3)
    imshow((single(squeeze(M.im(:,:,1))).*(0.5+mask/2))',[]);set(gca,'YDir','normal')
    title('combined')
    saveas(gcf,'Fig04 dose area mask.pdf')
    
    %check if nothing rotates in these areas
    figure(5);clf;
%     fullmask=uint16(~isnan(mask));
%     fullmask=uint16(repmat(fullmask,[1,1,size(M.im,3)]));
%     fullmask=fullmask.*M.im;

    for i=1:M.nFrame
       a=single(squeeze(M.im(:,:,i))).*mask;
       dosemin(i)=min(a(a>0));
    end
    plot(dosemin)
    ylabel('minimum image value of masked frames')
    xlabel('frame number')
    title('if this has major dips, then some object rotates into the dose mask')
     saveas(gcf,'Fig05 dose area purity check.pdf')
    
    
end


end

