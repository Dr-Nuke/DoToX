classdef lft
    
%% function collection file for the liquid film thickness computation
methods(Static) %evil function scope hack
    
function [c,rad,m]=FindCenters(F,reconcase)
    fprintf('case %d finding pin centers... \n',F.i)
    % finds the center point coordinates
    % F = F-struct
    % reconcase = XxXxZ Array of Z reconstructions of XxX pixel
    c=zeros(F.size(4),4,2); % center coordinates
    rad=zeros(F.size(4),4); % radius values found
    m=zeros(F.size(4),4);   % confidence measure
    
    for k=F.centeringrange%imrange(2,2) %iterate planes
        f.f_BoCount(k,20,20,5)
        %slice data
        %block(:,:,k)=imrotate(block(:,:,k),15,'crop');
        im=squeeze(reconcase(:,:,k));
        
        %create black/white-image (helps the algorhythm)
        imbw=~imbinarize(im,F.ImBinThresh(F.i));
        % find centers & radii
        try
            [c(k,:,:),rad(k,:),m(k,:)]=f_hugh_4_2(imbw,F.r_min,F.dr,F.S,F.o); %have this file in the same folder
            if ~all(rad(k,:))
                fprintf ('%4d zero! detected\n',k)
            end
        catch
            fprintf('%4d had a problem assigning the circle values\n',k)
        end
    end % k-loop
    fprintf('\n')
end

function [clinfit,cfit]=FitCenters(F)
    % finds a fit to the center coordinates to remove noise
    % first it removes "outliers" and performs a fit to the remaining
    % "good" points
    fprintf('case %d fitting centers... ',F.i)
    yplane=1:F.h;   % stack heigt or number of planes
    for i=1:F.geo.n_pins %iterate pins
        for j=1:2 % iterate x and y coordinate
            
            % find the median value of the center coordinates
            medicent(i,j)=median(squeeze(F.c(F.i,:,i,j)));
            
            % find those points that lie close to the median
            choice=abs(squeeze(F.c(F.i,:,i,j)-medicent(i,j)))<5;
            
            % apply to both yplane and center points
            ychoice=yplane(choice);
            cchoice=F.c(F.i,choice,i,j);
            
            %fit those values
            clinfit{i,j}=polyfit(ychoice,cchoice,1);
            
            % evaluate for the full stack
            cfit(:,i,j)=polyval(clinfit{i,j},yplane);
       end
    end
    fprintf('done.\n');
end

function checkCenters(F)
h=figure(18);clf
% check if the centers were found sufficiently good, and make the fit.
% plot everything
col=hsv(8);
for i=1:F.Ncases % case
    ax(i)=subplot(2,2,i);
    hold on
    for j=1:4 %pin
        p(2*j-1)=plot(squeeze(F.c(i,:,j,1)),'.','MarkerSize',5,'color',col(2*j-1,:),'DisplayName',sprintf('pin %d x',j));
        p(2*j)  =plot(squeeze(F.c(i,:,j,2)),'.','MarkerSize',5,'color',col(2*j,:),'DisplayName',sprintf('pin %d y',j));
        for k=1:2
            p(9)=plot(1:1024,squeeze(F.cfit(i,:,j,k)),'k','DisplayName','fit');
                            
            %plot([300,1020],[F.medicent(i,j,k),F.medicent(i,j,k)]+5,'k--')
            %plot([300,1020],[F.medicent(i,j,k),F.medicent(i,j,k)]-5,'k--')
        end
    end
    if or(i==1,i==3)
            ylabel('center point coordinate')
    end
        if or(i==2,i==3)
            xlabel('y-plane')
    end
    title(sprintf('rod center coordinates case %d',i))
    legend(p,'location','northwest')
    grid on
   
end
fname=sprintf('Fig18 pin center fits');
f.f_BoFig2PDF(h,fname)
end
  
function compareCenters(F)
    h=figure(20);clf;
    ax=gca;
    hold on
    col=hsv(3);
    for i=1:F.Ncases
        for j=1:F.geo.n_pins
            for k=1:2
                p(k)=plot(squeeze(F.cfit(i,:,j,k)),'color',col(k,:),...
                    'DisplayName',sprintf('c %d, p %d',i,j));
            
            end
        end
    end
    legend(p,'location','east')
    grid on
    xlabel('y-plane')
    ylabel('center coordinate')
    title('center coordinates across cases')
    fname=sprintf('Fig20 center comparison');
    f.f_BoFig2PDF(h,fname)
end

function [ProfEndPnt,ProfAngles]=FindProfilePaths(F)
    % creates the profile end points
    fprintf('case %d creating improfile coordinates... ',F.i)
    ProfEndPnt=zeros(F.size(4),F.geo.n_pins,F.pth.n_angles,2);
    a=1:F.geo.n_pins; % pin indicator
    b=circshift(a,-1); %shifted indicator for pin_n minus pin(n-1)
    
    %starting angles
    a_start=squeeze(rad2deg(atan2(F.cfit(F.i,:,b,2)-F.cfit(F.i,:,a,2),...  %y-difference
        F.cfit(F.i,:,b,1)-F.cfit(F.i,:,a,1)))... %x-difference
        -(F.pth.d_angle-90)/2); % add the more than 90deg section
    
    % the angle map
    ang=linspace(0,F.pth.d_angle,F.pth.n_angles);% generate incremental angles first
    
    ProfAngles= degtorad(repmat(a_start,[1,1,F.pth.n_angles])... % start angles
        +repmat(reshape(ang,[1,1,F.pth.n_angles]),F.h,4,1)); % plus incremental angles
    
    % x values of end points
    ProfEndPnt(:,:,:,1) = squeeze(repmat(F.cfit(F.i,:,:,1),[1,1,1,F.pth.n_angles]))... % start at center x
        +F.pth.r_path*cos(ProfAngles);     % and move r_path into Profangles direction
    ProfEndPnt(:,:,:,2) = squeeze(repmat(F.cfit(F.i,:,:,2),[1,1,1,F.pth.n_angles]))... % start at center x
        +F.pth.r_path*sin(ProfAngles);     % and move r_path
    
    fprintf('done.\n')
end

function [ProfEndPnt,ProfAngles]=FindProfilePathsM(F)
    % creates the profile end points for Michas tomos
    %fprintf('case %d creating improfile coordinates... ',F.i)
    ProfEndPnt=zeros(F.geo.n_pins,F.pth.n_angles,2);
    a=1:F.geo.n_pins; % pin indicator
    b=circshift(a,-1); %shifted indicator for pin_n minus pin(n-1)
    
    %starting angles
    a_start=squeeze(rad2deg(atan2(F.cen.cen(b,2)-F.cen.cen(a,2),...  %y-difference
        F.cen.cen(b,1)-F.cen.cen(a,1)))... %x-difference
        -(F.pth.d_angle-90)/2); % add the more than 90deg section
    
    
    
    % the angle map
    ang=linspace(0,F.pth.d_angle,F.pth.n_angles);% generate incremental angles first
    
    ProfAngles= degtorad(repmat(a_start,[1,F.pth.n_angles])... % start angles
        +repmat(reshape(ang,[1,F.pth.n_angles]),F.geo.n_pins,1)); % plus incremental angles
    
    % x values of end points
    

    ProfEndPnt(:,:,1) = squeeze(repmat(F.cen.cen(:,1),[1,F.pth.n_angles]))... % start at center x
        +F.pth.r_path*cos(ProfAngles);     % and move r_path into Profangles direction
    ProfEndPnt(:,:,2) = squeeze(repmat(F.cen.cen(:,2),[1,F.pth.n_angles]))... % start at center x
        +F.pth.r_path*sin(ProfAngles);     % and move r_path
    
   
end


function ImProfiles=FindProfiles(F,recon)
    ImProfiles=zeros(F.h,F.geo.n_pins,F.pth.n_angles,F.pth.r_path+1);

    % parfor doesnt like structs, so i copy the values needed
    h=F.h;
    n_pins=F.geo.n_pins;
    n_angles=F.pth.n_angles;
    cfit=squeeze(F.cfit(F.i,:,:,:));
    ProfEndPnt=squeeze(F.ProfEndPnt(F.i,:,:,:,:));
    r_path=F.pth.r_path+1;

    fprintf('case %d getting ImProfiles... \n',F.i)
    
    
    parfor i=1:h %planes
        f.f_BoCount(i,20,10,5)
         for j=1:n_pins 
            for k=1:n_angles 
                try
                    ImProfiles(i,j,k,:)=...
                        reshape(improfile(recon(:,:,i)',... % profile of the image
                        [cfit(i,j,1),ProfEndPnt(i,j,k,1)],... % from start point
                        [cfit(i,j,2),ProfEndPnt(i,j,k,2)],... % to end point
                        r_path),...       % with these many samples
                        [1,1,1,r_path]);    % and play with reshaping abit
                catch
                    fprintf('cought %d %d %d %d %d',F.i,i,j,k)
                end
            end
        end
    end
    
    fprintf('done.\n')
end

function ImProfiles=FindProfilesM(F,recon)
    ImProfiles=zeros(max(F.raw.nplane),F.geo.n_pins,F.pth.n_angles,F.pth.r_path+1);

    % parfor doesnt like structs, so i copy the values needed
    n_pins=F.geo.n_pins;
    n_angles=F.pth.n_angles;
    h=F.h;
    cfit=F.cen.cen;
    ProfEndPnt=F.Pa.ProfEndPnt;
    r_path=F.pth.r_path+1;
    fprintf('getting ImProfiles...')
    for pla=1:h %planes
        %disp(pla)
        
         for pin=1:n_pins 
            for ang=1:n_angles 
                try
                    ImProfiles(pla,pin,ang,:)=...
                        reshape(improfile(recon(:,:,pla)',... % profile of the image
                        [cfit(pin,1),ProfEndPnt(pin,ang,1)],... % from start point
                        [cfit(pin,2),ProfEndPnt(pin,ang,2)],... % to end point
                        r_path),...       % with these many samples
                        [1,1,1,r_path]);    % and play with reshaping abit
                catch
                    sprintf('cought  %d %d %d %d',pla,pin,ang)
                end
            end
        end
    end
    fprintf('done.\n')
end

function ProfileAnalyser(F,recon)
    
    planes=[1010,450,350,90];
    lpm=[0 3.45,2.85,2.3]
    pin=1;
    col=jet(F.pth.n_angles-6);
    h=figure(21);clf;
    ang=linspace(0,F.pth.d_angle,F.pth.n_angles);% generate incremental angles first
    angles=ang-90; % convert to [-90,90]
    for cas=1:4
        for planeind=1:length(planes)
            plane=planes(planeind);
            axind=(planeind-1)*F.Ncases+cas;
            ax(axind)=subplot(length(planes),F.Ncases,axind) ;           
            hold on
            for ii=1:(F.pth.n_angles-6) % the 2x3 outter angles are not "clean"
                i=ii+3;
                try
                    p(ii)=plot(squeeze(F.pth.Profs(cas,plane,pin,i,:)),...
                        'color',col(ii,:),'Displayname',sprintf('%d°',angles(i)));
                catch
                    fprintf('cas %d plane %d pin %d i %d\n',cas,plane,pin,i)
                    disp(size(F.pth.Profs))
                end
            end
            grid on
            if axind > 9
                xlabel('distance to pin center, [pixel]')
            end
            if ~all([1 4 7 10]-axind)
                ylabel('grey value')
            end
            if axind ==1
                legend([p(8:20:(F.pth.n_angles-6))],'location','northwest')
            end
            titlestring=sprintf('%.2f l/m, plane %d, pin %d\n',lpm(cas),plane,pin');
            if axind==2
                titlestring=sprintf('angular profile comaprison across cases\n %s',...
                    titlestring);
            end
            title(titlestring)
                
            ylim(0.0001*[-2,12])
        end
        linkaxes(ax)
    end
    fname=sprintf('Fig20 angular profile overview',cas, plane);
%    f.f_BoFig2PDF(h,fname)
    
    fid=211
    h=figure(fid);clf;
    for cas=1:4
        for planeind=1:length(planes)
            plane=planes(planeind);
            axind=(planeind-1)*F.Ncases+cas;
            ax(axind)=subplot(length(planes),F.Ncases,axind) ;           
            hold on

                try
                    
                    imagesc(squeeze(F.pth.Profs(cas,plane,pin,:,:)))
%                     p(ii)=plot(squeeze(F.pth.Profs(cas,plane,pin,i,:)),...
%                         'color',col(ii,:),'Displayname',sprintf('%d°',angles(i)));
                catch
                    fprintf('cas %d plane %d pin %d \n',cas,plane,pin)
                    disp(size(F.pth.Profs))
                end
            

            if axind > 9
                xlabel('distance to pin center, [pixel]')
            end
            if ~all([1 4 7 10]-axind)
                ylabel('angle')
            end
            
            titlestring=sprintf('%.2f l/m, plane %d, pin %d\n',lpm(cas),plane,pin');
            if axind==2
                titlestring=sprintf('angular profile comaprison across cases\n %s',...
                    titlestring);
            end
            title(titlestring)
                
            
        end
        linkaxes(ax)
    end
    fname=sprintf('Fig%d angular profile overview',fid);
%    f.f_BoFig2PDF(h,fname)
    
    
end

function moreplots(F,recon)
planes=[1010,450,350,90];
    lpm=[3.45,2.85,2.3]
    pin=1;
    col=jet(F.pth.n_angles-6);
    h=figure(22);clf;
    ang=linspace(0,F.pth.d_angle,F.pth.n_angles);% generate incremental angles first
    angles=ang-90; % convert to [-90,90]
    for cas=1:3
        for planeind=1:length(planes)
            plane=planes(planeind);
            axind=(planeind-1)*3+cas;
            ax(axind)=subplot(4,3,axind) ;           

                try
                    imshow(squeeze(recon(cas,:,:,plane))',0.0001*[-5,15])
                    set(gca,'YDir','normal')
                catch
                    fprintf('error case %d plane %d pin %d\n',cas,plane,pin)
                    disp(size(F.pth.Profs))
                end
            axis tight
            grid on



            titlestring=sprintf('%.2f l/m, plane %d, pin %d\n',lpm(cas),plane,pin');
            if axind==2
                titlestring=sprintf('tomogram overview\n %s',...
                    titlestring);
            end
            title(titlestring)

        end
    end
    fname=sprintf('Fig22 planes overview',cas, plane);
    f.f_BoFig2PDF(h,fname)
    
end

function bg= d2oBackground(F,recon)
    % calculates the d2o background of a reconstruction
    bg=zeros(F.Ncases,F.h,F.geo.n_pins); %BackGround
    
    r_bg=30; % radius up to which the background calculation shall be based on
    % calculate the volume weighted average and divide by the Area
    dr=repmat(reshape(0:r_bg,[1,1,1,1,r_bg+1]),[1,F.h,F.geo.n_pins,F.pth.n_angles,1]);
    bg=mean(sum(F.pth.Profs(F.i,:,:,:,1:(r_bg+1)).*dr,5)./sum(dr,5),4);
end

function CheckBg(F)
    % checks if the background is reasonable
    h=figure(23);clf;
    ax=gca();
    for i=1:F.Ncases
        subplot(2,2,i)
        hold on;
        for j=1:F.geo.n_pins;
            plot(squeeze(F.bgHeat(i,:,j)),'displayname',sprintf('c%d p%d',i,j));
        end;
        ylim(0.0001*[0,1])
        xlim([1,1024])
        grid on
        legend()
        xlabel('y-plane')
        ylabel('mean d2o area grey level')
        title(sprintf('case %d background',i))
    end
    fname=sprintf('Fig23 grayvalue background of heating area');
    f.f_BoFig2PDF(h,fname)
end


function lft=LFT(F)
    % finally calculating the LFT signal
    lft=sum(F.pth.Profs(:,:,:,:,F.lftSumRange),5);%... % sum up the grey value
end

function LFTMap(F)
    
    clim=[0,120];
    for cas=1:F.Ncases
       h(cas)=figure(23+cas); 
       for pin=1:F.geo.n_pins
           subplot(1,F.geo.n_pins,pin)
           imshow(squeeze(F.LFT(cas,:,pin,:)),[]);set(gca,'YDir','normal')
           
           if pin==1;
               ylabel('rod height')
           end
           if pin==F.geo.n_pins
               cb=colorbar();
               ylabel(cb, 'LFT [$\mu m$]','Interpreter','latex')
           end
           axis on
           set(gca,'ytick',[])
           set(gca,'xtick',[])
           xlabel('ang. pos.')
           title(sprintf('pin %d',pin))
           caxis(clim)
       end
    fname=sprintf('Fig%d LFT maps per case %d',cas+23,cas);
    f.f_BoFig2PDF(h(cas),fname)
    end
    
    for pin=1:2
        h(pin)=figure(26+pin);clf;
        for cas=2:F.Ncases
           subplot(1,F.Ncases-1,cas-1)
           imshow(squeeze(F.LFT(cas,:,pin,:)),[]);set(gca,'YDir','normal')
           if pin==1;
               ylabel('rod height')
           end
           if cas==F.Ncases
               cb=colorbar();
               ylabel(cb, 'LFT [$\mu m$]','Interpreter','latex')
           end
           axis on
           set(gca,'ytick',[])
           set(gca,'xtick',[])
           xlabel('ang. pos.')
           title(sprintf('%.2f l/min',F.flow(cas)))
           caxis(clim)
        end
        fname=sprintf('Fig%d LFT maps per pin %d',pin+26,pin);
        f.f_BoFig2PDF(h(cas),fname)
    end
    
    
end

function plotProfiles(F)
    
    figure(30);clf;
    ax=gca();
    hold on
    col=hsv(F.h);
    x=zeros(1,F.pth.r_path+1);
    
    cas=2;
    pin=1;
    ang=91;
    for pla=1:3:F.h
       plot3(1:F.pth.r_path+1,x+pla-1,squeeze(F.pth.Profs(cas,pla,pin,ang,:)))%-...
%           F.bgHeat(cas,pla,pin),'color',col(pla,:))
    end
    xlabel('pin center distance [pixel]')
    ylabel('channel height [pixel]')
    zlabel('grey value')
    xlim([0,F.pth.r_path+1])
    ylim([0,F.h])
    zlim(0.001*[-0.5,1.5])
    grid on
    title('grey value profiles ')
    view([37,50])
    
    
    figure(31);clf;
    % doesnt look good
    ax=gca();
    x=0:F.pth.r_path;
    y=1:F.h;
    [xx,yy]=meshgrid(x,y);
    surf(xx,yy,squeeze(F.pth.Profs(cas,:,pin,ang,:)))%-...
         %  repmat(reshape(F.bgHeat(cas,:,pin),[F.h,1]),[1,F.pth.r_path+1]))
    xlabel('pin center distance [pixel]')
    ylabel('channel height [pixel]')
    zlabel('grey value')
    xlim([0,F.pth.r_path+1])
    ylim([0,F.h])
    zlim(0.001*[-0.5,1.5])
    grid on
    title('grey value profiles ')
    
    
    % horizontal profiles
    for pin=[1,2]
    h=figure(310+pin);clf; 
    set(gcf,'pos',[10 10 600 900])
    cas=2;

    ax=gca(); hold on;
    col=hsv(F.Ncases);
    col=[1 1 1;
        1 0 0;
        0 1 0;
        0 0 1];
    dist=1;
    Nplots=13;
    heights=linspace(0,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix(heights); % calculate round pixel planes
    angles=linspace(-90,90,181);

    for cas=2:F.Ncases
        for i=1:length(planes)
            pla=planes(i);
            if i==1;
                
                p(cas)=plot(angles,squeeze(mean(F.LFT(cas,10:19,pin,:),2))+dist*0,...
                'color',col(cas,:),'DisplayName',...
                sprintf('case %d, %.1f l/min ',cas,F.flow(cas)));
            line([-90,90],0*[pla,pla],'color',0.5*[1,1,1])
            else
            
            p(cas)=plot(angles,squeeze(mean(F.LFT(cas,pla:pla+9,pin,:),2))+dist*pla,...
                'color',col(cas,:),'DisplayName',...
                sprintf('case %d, %.1f l/min',cas,F.flow(cas)));
            line([-90,90],dist*[pla,pla],'color',0.5*[1,1,1])
            end
        end
    end
    
    for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    % add spacer indicator
    spacery=([165,290]);%lft.pix2hcm([165,290]);
    x = 90*[-1,1,1,-1];
    y = 0.015*[1,1,2,2];
    y=dist*[spacery(1),spacery(1),spacery(2),spacery(2)]
    fill(x, y, 'k','facealpha',0.2,'edgecolor','none')
    text(60,dist*230,'spacer','color','k')
    
    % add vanes indicator
    vanesy=([290,331]);
    x = 90*[-1,1,1,-1];
    y=dist*[vanesy(1),vanesy(1),vanesy(2),vanesy(2)]
    fill(x, y, [0,0,1],'facealpha',0.2,'edgecolor','none')
    text(60,310*dist,'vanes','color','b')
    
    grid on
    yticks(planes(1:2:end)*dist)
    yticklabels(ylab(1:2:end))
    %     yticklabels(ylab)
    xlim([-90,90])
    xlabel('angular position [deg]')
    xticks([-90,-45,0,45,90])
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    
    ylabel('channel length [cm]')
    ylim([0,1080*dist])
    title(sprintf('horizontal LFT profiles pin %d',pin))
    legend(p(2:F.Ncases),'location','southoutside')

    fname=sprintf('Fig33 LFT profiles pin %d.png',pin);    
    set(gcf,'color','w')
    set(gcf,'position',[-1145 -9 326 901])
    saveas(gcf,fname)
    %lft.f_BoFig2PDF(h,fname)
    end
    
end

function plotMinMaxProfiles(F)
    
    % horizontal profiles
    for pin=1:2
    h=figure(35+pin);clf; 
    ax=gca();
    hold on
    cas=1;
   
    ang=4:178;
        col=[1 1 1;
        1 0 0;
        0 1 0;
        0 0 1];
    col2=1-col;
    Nplots=14;
    heights=linspace(-1,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix(heights); % calculate round pixel planes


    for cas=2:F.Ncases
        p(cas,1)=plot(squeeze(max(F.LFT(cas,:,pin,ang),[],4)),'color',col(cas,:),'displayname','min');
        p(cas,2)=plot(squeeze(mean(F.LFT(cas,:,pin,ang),4)),'color',col(cas,:)+0.4*col2(cas,:),'displayname','max');
        p(cas,3)=plot(squeeze(min(F.LFT(cas,:,pin,ang),[],4)),'color',col(cas,:)+0.7*col2(cas,:),'displayname','mean');
    end
    pp=plot([-200,-100],[0 0],'r');
    %text(300,0.008,sprintf('pin %d'),pin)
    angles=linspace(-90,90,181);
    grid on;
    xlabel('channel height')
    ylabel('LFT, [um]')
    xlim([1,F.h])
    yli=[-20,150];
    ylim(yli)
    xticks(planes)
    for tick=1:length(heights)
        xt{tick}=num2str(heights(tick));
    end
    xticklabels(xt)
    
    legend([p(2:4,1)',p(2,2),p(2,3)],{'3.5 l/min','2.9 l/min','2.3 l/min','max','mean','min'},'location','northeast')
    view([90 -90])

    line([0,F.h],[0,0],'color',[0 0 0])
    fname=sprintf('Fig%d minmax plot pin %d.png',35+pin,pin);
    saveas(gcf,fname)
    %lft.f_BoFig2PDF(h,fname)
    end
    
    
end

function Rayplot(F,recon)
    % plots the reconstruction with rays
    h=figure(34);clf;
    cas=1;
    pin=1;
    pla=667;
    imshow(squeeze(sum(recon(cas,:,:,pla:pla+9),4))',[]);
    set(gca,'YDir','normal')
    hold on
    col=jet(F.pth.n_angles);
    for pin=1:F.geo.n_pins
    for ang=1:10:F.pth.n_angles
        plot([F.cfit(cas,pla,pin,1),F.ProfEndPnt(cas,pla,pin,ang,1)],...
             [F.cfit(cas,pla,pin,2),F.ProfEndPnt(cas,pla,pin,ang,2)],...
             'color',col(ang,:))
    end
    end

    

end

function circleplot(F,recon)
    h=figure(35);clf
    cas=2;
    pla=481;
    imshow(squeeze(recon(cas,:,:,pla))',[]);
    
    set(gca,'YDir','normal');
    hold on
    
    %     fname=sprintf('Fig35 reconstruction');
    %     lft.f_BoFig2PDF(h,fname)
    for pin=1:F.geo.n_pins
        viscircles(squeeze(F.c(cas,pla,pin,:))',squeeze(F.rad(cas,pla,pin)),...
            'LineStyle','--','EnhanceVisibility',0)
        
    end
    
    
    
    scatter(squeeze(F.c(cas,pla,:,1)),squeeze(F.c(cas,pla,:,2)),'r+','linewidth',1)
    scatter(squeeze(F.cfit(cas,pla,:,1)),squeeze(F.cfit(cas,pla,:,2)),'gx','linewidth',1)
    %     fname=sprintf('Fig35 circle plot');
    %     lft.f_BoFig2PDF(h,fname)
    
    col=cool(F.pth.n_angles)
    for pin=1:F.geo.n_pins
        for ang=1:10:F.pth.n_angles
            plot([F.cfit(cas,pla,pin,1),F.ProfEndPnt(cas,pla,pin,ang,1)],...
                [F.cfit(cas,pla,pin,2),F.ProfEndPnt(cas,pla,pin,ang,2)],...
                'color',col(ang,:))
        end
    end

%     fname=sprintf('Fig35 fan plot');
%     lft.f_BoFig2PDF(h,fname)
    
    h=figure(36);clf;
    ax=gca();
    hold on
    col=jet(F.pth.n_angles);
    for ang=1:3:F.pth.n_angles
       p(ang)=plot(squeeze(F.pth.Profs(cas,pla,pin,ang,:)),'color',col(ang,:),...
       'displayname',sprintf('%d°',ang-91));
    end
    xlabel('pin center distance [pixel]')
    ylabel('channel height [pixel]')
    legend(p([1:30:end]), 'location','northwest');
    grid on
    fname=sprintf('Fig35 fan plot');
    lft.f_BoFig2PDF(h,fname)
end

function cm=pix2hcm(y)
    % converts the y plane into a cm value
    cm=4.5/(386-11)*y;   % data taken from 4.5cm wedge gauge and the according
    % half width pixels
end

function cm=pix2hcm0(y)
    % converts the y plane into a cm value, with the spacer beginning as
    % zero
    cm=4.5/(386-11)*(y-163);   % data taken from 4.5cm wedge gauge and the according
    % half width pixels
end

function y=hcm2pix(cm)
   % converts a cm to pixels
   y=round((386-11)/4.5*cm);
end

function y=hcm2pix0(cm)
   % converts a cm to pixels, with zero at spacer beginning
   y=round((386-11)/4.5*cm+163);
end

function dh=pix2dh(pix)
   % converts a pixel to hydraulic diameter
   dh=4.5/(386-11)*(pix-163)
end

function pix=dh2pix(dh)
   % converts a cm to pixels
   pix=round();
end

function []= f_BoFig2PDF(h,fname)
    % prints a figure to .fig,.pdf and .png
    % removed resizing
    set(h,'PaperOrientation','portrait');
    set(h,'PaperPosition', [1 1 19 28]);
    savefig(h,fname)
    print(h, '-dpdf',fname);
    set(h,'PaperOrientation','portrait');
    print(h, '-dpng',fname);
end

function reconFilter(F,recon)
    sigma = 5;
    sz = 30;    % length of gaussFilter vector
    x = linspace(-sz / 2, sz / 2, sz);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    gaussFilter = ones(size(x))/sz
    figure(405);clf;
    plot(x,gaussFilter)
    grid on

end

function [x]=coords(F)

numxmax1=(F.recsize-1)/2; % maximum pixel number in each direction
numx1=linspace(-numxmax1,numxmax1,F.recsize); % number of each pixel
coxmax1=F.recres*(F.recsize-1)/2; % coordinate of the furthest pixel
x=linspace(-coxmax1,coxmax1,F.recsize); % coordinate list
end

function mm=plane2mm(plane)
    % take plane number and get back mm
    mm=(plane-292)*0.127;
end

function plane=mm2plane(mm)
    plane=mm/0.127+292
end

function dh=mm2dh(mm)
    dh=mm/(8.473);
end

function dh=plane2dh(plane)
    dh=lft.mm2dh(lft.plane2mm(plane));
end

function plane=dh2plane(dh)
    plane=dh*(8.473)/0.127+292;
end

function rplane=plane2rplane(plane)
    rplane=(plane-10)/16+1;
end

function plane=rplane2plane(rplane)
    plane=(rplane+1)*16+10;
end


end %static
end %class
















