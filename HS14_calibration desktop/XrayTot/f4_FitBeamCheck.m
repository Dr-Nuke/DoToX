function []= f4_FitBeamCheck(M)
% checks if the fitbema worked properly

if M.proofs==1
mask=M.mask;
maskmesh=size(mask);
[x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2));
xx=x;
yy=y;
xx(mask==0)=[];
yy(mask==0)=[];

    % check if the fits are ok

    h=figure(6);clf;
    %check some random frames
    imax=10;

    ii=round(M.nFrame(M.i)*rand(1,imax));
    for i = 1:length(ii)

        im=squeeze(M.im(:,:,ii(i)));
        im(M.mask==0)=[];
        rf =10; % reduction factor

        plot(M.f{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
        hold on
        %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
        title(sprintf('%02d nun-uniformity data \n according to %d random frames',M.i,imax ))
        xlabel('x')
        ylabel('y')
        zlabel('counts')
        pause(0.5)
    end
    
    fname=sprintf(sprintf('Fig06 %02d Beam fit',M.i));
    f_BoFig2PDF(h,fname) 
    
    %check if the corrected values are ok
    
    h=figure(7);clf;
    subplot(2,3,4:5);
    hold on
    ndf=size(M.dose,2); % number of dose fields
    
    doseprofile=zeros(M.nFrame(M.i),ndf+1);
    oldprofile=zeros(M.nFrame(M.i),ndf+1);
    dosepix=zeros(1,ndf);

    for i=1:ndf %iterate dose fields
        % mean of each field for each frame, for both old (im) and
        % corrected (imc) data
        doseprofile(:,i)=squeeze(mean(mean((M.imc(M.dose(1,i):M.dose(2,i),...
            M.dose(3,i):M.dose(4,i),:)))));
        oldprofile(:,i)=squeeze(mean(mean((M.im(M.dose(1,i):M.dose(2,i),...
            M.dose(3,i):M.dose(4,i),:)))));
        
        % plot
        plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))
        
        % the number of dose mask pixels
        dosepix(i)=(M.dose(2,i)-M.dose(1,i)+1)*(M.dose(4,i)-M.dose(3,i)+1);
        
        % the total dose, should be constantly equal 1
        doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
        oldprofile(:,ndf+1)=oldprofile(:,ndf+1)+dosepix(i)*oldprofile(:,i);
    end
    %total dose manually done
    doseprofile(:,ndf+1)=doseprofile(:,ndf+1)/sum(dosepix);
    oldprofile(:,ndf+1)=oldprofile(:,ndf+1)/sum(dosepix);
    plot(doseprofile(:,ndf+1),'DisplayName','combined')

    title(sprintf('relative doses of the %d regions \n after dose correction',ndf))
    xlabel('frame Nr')
    ylabel('relative dose')
    grid on

    legend()
  
    
    
    subplot(2,3,1:2);
    plot(oldprofile(:,ndf+1))
    ylabel('mean pixel value')
    title(sprintf('%02d mean per frame dose of raw footage',M.i))
    grid on
    
    
    
    subplot(2,3,[3,6]);
    imshow(squeeze(M.im(:,:,1))',[]);set(gca,'YDir','normal')
    hold on
    for i=1:ndf
        r=rectangle('Position',[M.dose(1,i),...
                              M.dose(3,i),...
                              M.dose(2,i)-M.dose(1,i),...
                              M.dose(4,i)-M.dose(3,i)],...
                              'Facecolor','r');
    end
    title('dose areas')
    fname=sprintf('Fig07 %02d dose area checks',M.i);
    f_BoFig2PDF(h,fname) 
    M.d.doseprofile(M.i,:)=oldprofile(:,ndf+1);
end

    
end

