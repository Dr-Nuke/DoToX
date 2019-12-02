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

    figure(6);clf;
    %check some random frames
    imax=10;

    ii=round(M.nFrame*rand(1,imax));
    for i = 1:length(ii)

        im=squeeze(M.im(:,:,ii(i)));
        im(M.mask==0)=[];
        rf =10; % reduction factor

        plot(M.f{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
        hold on
        %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
        title(sprintf('nun-uniformity data \n according to %d random frames',imax ))
        xlabel('x')
        ylabel('y')
        zlabel('counts')
        pause(0.5)
    end
    saveas(gcf,'Fig06 Beam fit.pdf')
    
    %check if the corrected values are ok
    
    figure(7);clf;
    subplot(2,1,1);
    hold on
    ndf=size(M.dose,2); % number of dose fields
    
    doseprofile=zeros(M.nFrame,ndf+1);

for i=1:ndf
    doseprofile(:,i)=squeeze(mean(mean((M.imc(M.dose(1,i):M.dose(2,i),...
        M.dose(3,i):M.dose(4,i),:)))));
    plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))
    dosepix(i)=(M.dose(2,i)-M.dose(1,i)+1)*(M.dose(4,i)-M.dose(3,i)+1);
    doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
end
    %total dose manually done
    doseprofile(:,ndf+1)=doseprofile(:,ndf+1)/sum(dosepix);
    plot(doseprofile(:,ndf+1),'DisplayName','combined')

    title(sprintf('relative doses of the %d regions \n after dose correction',ndf))
    xlabel('frame Nr')
    ylabel('relative dose')
    zlabel('counts')
    legend()
    
    
    subplot(2,1,2);
    imshow(squeeze(M.im(:,:,1))',[]);set(gca,'YDir','normal')
    hold on
    for i=1:ndf
        r=rectangle('Position',[M.dose(1,i),...
                              M.dose(3,i),...
                              M.dose(2,i)-M.dose(1,i),...
                              M.dose(4,i)-M.dose(3,i)],...
                              'Facecolor','r');
    end
    
    
    saveas(gcf,'Fig07 corrected doses.pdf')
end

% remove data not needed anymore






end

