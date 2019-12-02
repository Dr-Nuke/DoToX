function M = f4_centering(M)
% does the centering
rpx=M.rpx;

% sinograms
M.sino=M.im(rpx(1):rpx(2),rpx(3):rpx(4),:);
M.sinol=M.sino;
M.sinosize=size(M.sino);

% temporary!!!
%M=rmfield(M,'sino')% sinograms



%find the 180° image of 1st image
i180=ceil(M.sinosize(3)/2);
im180=squeeze(M.sino(:,:,i180));
M.shift=zeros(1,M.sinosize(2));
M.angles=linspace(0,360,M.sinosize(3)+1);
M.angles(end)=[];

for i = 1:M.sinosize(2) %classic center shift calculation
%     im=squeeze(M.sino(:,i,:));
%     [c,lags]=xcov(im(:,i),flipud(im180(:,i)),'coeff');
%     if any(isnan(c)) %debug
%         disp(sprintf('%4d contains NaN',i))
%         %ind_nan(k)=-1;
%         continue
% %     else
% %         disp(sprintf('%4d does not contain NaN',i))
% %         %ind_nan(k)=1;
%     end
%     [~,ind]=max(c);
%     M.shift(i)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
%         (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
%     
%     iml=squeeze(M.sinol(:,i,:));
%     M.sinoc{i}=f3_center(iml,M.shift(i)); %centering
%     M.sinoc{i}([1,end],:)=[];               % cropping
    try
        imc=circshift(squeeze(M.sinol(:,i,:)),60,2);
        recon=iradon(imc,M.angles,'spline','Hann',1,M.recsize);
        disp([i,size(recon)])
        M.recon(:,:,i)=recon; %hardcode after-recon-cropping
    catch
        disp(sprintf('%4d didnt reconstruct',i))
    end
    
    
end










if M.proofs==1
    
    % check: moving through the frames
yplane=1;
for i=1:floor(M.sinosize(3)/2)
    i1802=i+i180;
    im=squeeze(M.sino(:,:,i));
    im180=squeeze(M.sino(:,:,i1802));
    [c,lags]=xcov(im(:,yplane),flipud(im180(:,yplane)),'coeff');
    if any(isnan(c)) %debug
        disp(sprintf('%4d contains NaN',i))
        %ind_nan(k)=-1;
        continue
%     else
%         disp(sprintf('%4d does not contain NaN',i))
%         %ind_nan(k)=1;
    end
    [~,ind]=max(c);
    shift(i)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
end
    
    % combined sinogram study plot
    h=figure(8);clf;
    subplot(3,4,[1:3])
    a=squeeze(M.sino(:,50,:));
    imshow(a,[])
    title('A sinogram')
    colorbar
    
    subplot(3,4,4)
    hist(a(:),100)
    title('sinogram')
    grid on
    
    subplot(3,4,[5:7])
    b=squeeze(M.sinol(:,50,:));
    imshow(b,[])
    title('A log sinogram')
    colorbar

    subplot(3,4,8)
    hist(b(:),100)
    title('sinogram')
    grid on
    
    subplot(3,4,[9:11])
    c=b;
    c(b<0)=max(b(:));
    imshow(c,[])
    title('A log sinogram. bright dots represent <0 values')
    colorbar
    h=gcf;

    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'Fig08 sinogram studies.pdf');
    
    
   % plot of a sinogram with lines at 1 & 180°
   figure(9);clf;
   imshow(squeeze(M.sino(:,yplane,:))',[]); set(gca,'YDir','normal')
   line([1,M.sinosize(1)],[1,1])
   line([1,M.sinosize(1)],[i180,i180])
   title('two lines 180° apart')
   saveas(gcf,'Fig09 180 deg frames.pdf')
   
   % plot of the 180° study
   figure(10);clf;
  
   plot(squeeze(M.sino(:,yplane,1)),'k','Displayname','reference')
   hold on
   
   ii=linspace(-2,2,5);
   imax=length(ii);
   col=hsv(imax);
   for i = 1:imax
       plot(flipud(squeeze(M.sino(:,yplane,i180+ii(i)))),'color',col(i,:),...
           'Displayname',sprintf('%d',ii(i)))
   end
   legend()
   grid on
   xlabel('x position')
   ylabel('detector signal')
   title(sprintf('study on the correct 180° frame: manual deviation at plane %d',yplane))
   saveas(gcf,'Fig10 180 deg study.pdf')
   
   figure(11)
   plot(M.shift)
   xlabel('pixel plane (in flow direction)')
   ylabel('off-centerdness [pixel]')
   title('centering shift vs plane')
   grid on
   saveas(gcf,'Fig11 center shift .pdf')
   
   figure(12)
   plot(shift)
   xlabel('frame pair (frame i and frame i+180°)')
   ylabel('off-centerdness [pixel]')
   title('centering shift vs frame pair')
   grid on
   saveas(gcf,'Fig12 center shift .pdf')
   
   figure(13);clf;
   cropim=squeeze(M.im(:,:,1));
   mask=zeros(size(cropim));
   mask(rpx(1):rpx(2),rpx(3):rpx(4))=1;
   cropim=cropim+mask;
   imshow(cropim',[]);set(gca,'YDir','normal')
   title('reconstruction mask')
   saveas(gcf,'Fig13 reconstruction mask .pdf')
    
end



end

