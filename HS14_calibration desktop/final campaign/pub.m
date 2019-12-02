classdef pub
    
%% function collection file for the liquid film thickness computation
methods(Static) %evil function scope hack
    
function quadtomo(F,recon)
    h=figure(200);clf;
    cas=2;
    height=[-1,4,12,12];
    planes=lft.hcm2pix0(F.dh*height);
    xsh=[0,1,0,1].*size(recon,2); %shifts in x and y direction
    ysh=[1,1,0,0].*size(recon,2);
    
    %merge tomos
    for i=1:4
        pla=planes(i);
        im(:,:,i)=squeeze(mean(recon(cas,:,:,pla:pla+9),4));
    end
        im2=cat(2,im(:,:,4),im(:,:,3));
        im3=cat(2,im(:,:,2),im(:,:,1));
    im4=cat(2,im2,im3);
    imshow(im4',[]);set(gca,'YDir','normal');
    
    % add circles
    xshift=size(recon,2);
    for pin=1:F.geo.n_pins
        viscircles(squeeze(F.c(cas,pla,pin,:))'+[xshift,0],...
            squeeze(F.rad(cas,pla,pin)),...
            'LineStyle','--','EnhanceVisibility',0)
    end
    hold on
    
   fs=18
   col=spring(F.pth.n_angles);
   textshift=[-20,20,20,-20;
               -20 -10 10 20];
   xshift=0 ;
   yshift= size(recon,2)
   for pin=1:F.geo.n_pins
        for ang=1:18:F.pth.n_angles
            plot([F.cfit(cas,pla,pin,1),F.ProfEndPnt(cas,pla,pin,ang,1)]+xshift,...
                [F.cfit(cas,pla,pin,2),F.ProfEndPnt(cas,pla,pin,ang,2)],...
                'color',col(ang,:));
        end
        % pin labels
        text(F.cfit(cas,pla,pin,1)+textshift(1,pin)+xshift,...
            F.cfit(cas,pla,pin,2)+textshift(2,pin),...
            sprintf('%d',pin),'color',[0.5 0.5 1],'fontsize',fs)
    end
    % add textes
    % plane labels
    
    text(30,280+3*yshift,'$-d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30,280+2*yshift,'$4d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30,280+yshift,'$12d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30,280,'$12d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    
    
    % angle labesl
    text(F.ProfEndPnt(cas,pla,1,1,1)+xshift+4,F.ProfEndPnt(cas,pla,1,1,2)-5,...
        '$0$ ','Interpreter','latex','color',[1 0 1],'fontsize',fs-2)
    text(F.ProfEndPnt(cas,pla,1,end,1)+xshift-12,F.ProfEndPnt(cas,pla,1,end,2),...
        '$\pi$ ','Interpreter','latex','color',[1 1 0],'fontsize',fs-2)
    

    h.InvertHardcopy = 'off';
    fname=sprintf('Fig200 pub horizontal profiles');
    
    set(h,'color','w');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'filename','-dpdf','-r0')

    
    
%     u = h.Units;
%     h.Units = 'inches'
    export_fig test2.png
    h.PaperPositionMode = 'auto'
    h_pos = fig.PaperPosition;
    h.PaperSize = [h_pos(3) h_pos(4)];
    savefig(h,fname)
    
    print(h, '-dpdf',fname);

    print(h, '-dpng',fname);

end 

function quadtomo2x2(F,recon)
    fid=201
    h=figure(fid);clf;
    cas=2;
    height=[-1,4,12,12];
    planes=lft.hcm2pix0(F.dh*height);
    xsh=[0,1,0,1].*size(recon,2); %shifts in x and y direction
    ysh=[1,1,0,0].*size(recon,2);
    
    %merge tomos
    for i=1:4
        pla=planes(i);
        im(:,:,i)=squeeze(mean(recon(cas,:,:,pla:pla+9),4));
    end
        im2=cat(1,im(:,:,3),im(:,:,4));
        im3=cat(1,im(:,:,1),im(:,:,2));
    im4=cat(2,im2,im3);
    imshow(im4',[]);set(gca,'YDir','normal');
    
    % add circles
    xshift=size(recon,2);
    for pin=1:F.geo.n_pins
        viscircles(squeeze(F.c(cas,pla,pin,:))'+[xshift,0],...
            squeeze(F.rad(cas,pla,pin)),...
            'LineStyle','--','EnhanceVisibility',0)
    end
    hold on
    
   fs=18
   col=spring(F.pth.n_angles);
   textshift=[-20,20,20,-20;
               -20 -10 10 20];
   xshift=size(recon,2) ;
   yshift= size(recon,2)
   %add beams
   for pin=1:F.geo.n_pins
        for ang=1:18:F.pth.n_angles
            plot([F.cfit(cas,pla,pin,1),F.ProfEndPnt(cas,pla,pin,ang,1)]+xshift,...
                [F.cfit(cas,pla,pin,2),F.ProfEndPnt(cas,pla,pin,ang,2)],...
                'color',col(ang,:));
        end
        % pin labels
        text(F.cfit(cas,pla,pin,1)+textshift(1,pin)+xshift,...
            F.cfit(cas,pla,pin,2)+textshift(2,pin),...
            sprintf('%d',pin),'color',[0.5 0.5 1],'fontsize',fs)
    end
    % add textes
    % plane labels
    
    text(30,280+1*yshift,'$-d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30+xshift,280+1*yshift,'$4d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30,280,'$12d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    text(30+xshift,280,'$12d_h$ ','Interpreter','latex','color',[1 1 1],'fontsize',fs)
    
    
    % angle labesl
    text(F.ProfEndPnt(cas,pla,1,1,1)+xshift+4,F.ProfEndPnt(cas,pla,1,1,2)-5,...
        '$0$ ','Interpreter','latex','color',[1 0 1],'fontsize',fs-2)
    text(F.ProfEndPnt(cas,pla,1,end,1)+xshift-12,F.ProfEndPnt(cas,pla,1,end,2),...
        '$\pi$ ','Interpreter','latex','color',[1 1 0],'fontsize',fs-2)
    

    h.InvertHardcopy = 'off';
    fname=sprintf('Fig%d pub horizontal profiles',fid);
    
    set(h,'color','w');
    saveas(h,strcat(fname,'manual4x4.png'))
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'filename','-dpdf','-r0')

    
    
%     u = h.Units;
%     h.Units = 'inches'
    export_fig test2.png
    h.PaperPositionMode = 'auto'
    h_pos = fig.PaperPosition;
    h.PaperSize = [h_pos(3) h_pos(4)];
    savefig(h,fname)
    
    print(h, '-dpdf',fname);

    print(h, '-dpng',fname);

end 

function horzprof(F)
 % horizontal profiles
    for pin=[1,2]
    h=figure(203);clf; 
    set(gcf,'pos',[10 10 600 900])
  
    col=hsv(3);
    dist=0.0001;
    Nplots=15;
    heights=linspace(-2,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix0(heights*F.dh); % calculate round pixel planes
    angles=linspace(-90,90,181);
    
    
        ax(pin)=gca();
        hold on
    for cas=1:F.Ncases
        for i=1:length(planes)
            pla=planes(i);
%             if i==1;
%                 
%                 p(cas)=plot(angles,squeeze(mean(F.LFT(cas,10:19,pin,:),2))+dist*0,...
%                 'color',col(cas,:),'DisplayName',...
%                 sprintf('case %d, %.1f l/min ',cas,F.flow(cas)));
%             line([-90,90],0*[pla,pla],'color',0.5*[1,1,1])
%             else
            
            p(cas)=plot(angles,squeeze(mean(F.LFT(cas,pla:pla+9,pin,:),2))+dist*pla,...
                'color',col(cas,:),'DisplayName',...
                sprintf('case %d, %.1f l/min',cas,F.flow(cas)));
            line([-90,90],dist*[pla,pla],'color',0.5*[1,1,1])
%             end
        end
    end
    
    for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    % add spacer indicator
    spacery=([163,290]);%lft.pix2hcm([165,290]);
    x = 90*[-1,1,1,-1];
    y=dist*[spacery(1),spacery(1),spacery(2),spacery(2)];
    fill(x, y, 'k','facealpha',0.3,'edgecolor','none')
    text(30,0.022,'spacer','color','k')
    
    % add vanes indicator
    vanesy=([290,331]);
    x = 90*[-1,1,1,-1];
    y=dist*[vanesy(1),vanesy(1),vanesy(2),vanesy(2)];
    fill(x, y, [0,0,1],'facealpha',0.2,'edgecolor','none')
    text(30,0.031,'vanes','color','b')
    ax(pin).XGrid = 'on';
    if 1%pin ==1
        yticks(planes(1:2:end)*dist)
        yticklabels(ylab(1:2:end))
        ylabel('channel length [$d_h$]','Interpreter','latex')
    else
        yticks([])
    end
    xlim([-90,90])
    xlabel('angular position [rad]')
    xticks([-90,-45,0,45,90])
    xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
    ylim([0,0.108])
    title(sprintf('Pin %d',pin))
    
    hL=legend(p(1:3),'location','southoutside')
    
    fname=sprintf('Fig203 pub horizontal profiles pin %d',pin);
    set(h,'position',[0 0 305 820])
    
%     set(h,'color','w');
%     set(h,'Units','Inches');
%     pos = get(h,'Position');
%     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
     print(h,fname,'-dpdf','-r0')
     print(h,fname,'-dpng','-r0')
%     lft.f_BoFig2PDF(h,fname)
    end
    
    
    % lets try the 360°
    h=figure(204);clf; 
    set(gcf,'pos',[10 10 600 900])
  
    col=hsv(3);
    dist=0.0001;
    Nplots=15;
    heights=linspace(-2,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix0(heights*F.dh); % calculate round pixel planes

    angb=linspace(-50,50,101); %angle basis
    angi=angb+90;             % angle index for array
   
    ax=gca();    
    hold on
    for cas=1:F.Ncases
        for i=1:length(planes)
            pla=planes(i);
            for pin=[1,2,3,4]
                %disp([cas,i,pla,pin])
                angx=(-5:95)+90*(pin-1); % angle in x vector forplot
                p(cas)=plot(angx,squeeze(mean(F.LFT(cas,pla:pla+9,pin,angi),2))+dist*pla,...
                    'color',col(cas,:),'DisplayName',...
                    sprintf('case %d, %.1f l/min',cas,F.flow(cas)));
                line([0,360],dist*[pla,pla],'color',0.5*[1,1,1])
            end
%             end
        end
    end
    
    for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    
    for pin=[1,2,3,4]
        text( 30+ (pin-1)*90,0.107,sprintf('Pin %d',pin),'color','k')
    end
    
    yticks(planes(1:2:end)*dist)
    yticklabels(ylab(1:2:end))
    ylabel('channel length [$d_h$]','Interpreter','latex')

    ylim([0,0.108])
    
    % add spacer indicator
    spacery=([163,290]);%lft.pix2hcm([165,290]);
    x = 360*[0,1,1,0];
    y=dist*[spacery(1),spacery(1),spacery(2),spacery(2)];
    fill(x, y, 'k','facealpha',0.3,'edgecolor','none')
    text(220,0.022,'spacer','color','k')
    
    % add vanes indicator
    vanesy=([290,331]);

    y=dist*[vanesy(1),vanesy(1),vanesy(2),vanesy(2)];
    fill(x, y, [0,0,1],'facealpha',0.2,'edgecolor','none')
    text(220,0.032,'vanes','color','b')
    ax.XGrid = 'on';

    xlim([0,360])
    xlabel('angular position [rad]')
    xticks([0,90,180,270,360])
    xticklabels({'0','\pi/2','\pi','3/2\pi','2\pi'})

    title(sprintf('merged LFTP profiles'))
    legend(p(1:3),'location','southoutside')
    
    

fname=sprintf('Fig204 pub horizontal profiles 360');
        set(h,'color','w');
    saveas(h,strcat(fname,'manual.png'))
    lft.f_BoFig2PDF(h,fname)
    
    


end

function dhplot(F,d)
    fid=205;h=figure(fid);clf;
    ax=gca();
    imshow(squeeze(d(180:480,:,250))');set(gca,'YDir','normal')
    ylim([5,1020])
    fs=18;
    Nplots=15;
    heights=linspace(-2,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix0(heights*F.dh); % calculate round pixel planes
        for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    % add spacer indicator


    %text(30,0.022,'spacer','color','k')
    yticks()
  
    set(gca,'YTick',planes(1:2:end), 'YTickLabel', ylab(1:2:end), 'fontsize', fs);
    axis on
    xticks([])
    text(103,230,'spacer','fontsize',fs,'color',[1 1 1])
    annotation('textarrow',[90,95]/200,[420,370]/1020,'String','vanes',...
        'fontsize',fs)
    angau=sprintf('angular\n gauge');
    annotation('textarrow',[100,130]/200,[150,120]/1020,'string',...
        angau,'fontsize',fs-2)
    ylabel(sprintf('channl height [$d_h$]'),'Interpreter','latex','fontsize', fs)
    

    fname=sprintf('Fig%d projection',fid);
    
    lft.f_BoFig2PDF(h,fname)
end

function axprofile(F)
    fid=206;h=figure(fid);clf;
    ax=gca();
    hold on
    ang=45:135;
    col=hsv(3);
    colmod=[1,0.5,1,0.5]
    mark={}
    for cas=1:F.Ncases
        for pin=3:4%F.geo.n_pins
            p(cas,pin)=plot(mean(F.LFT(cas,:,pin,ang),4),'color',...
                colmod(pin)*col(cas,:),...
                'displayname', sprintf('case %d pin %d',cas,pin));
            
        end
    end
    legend()
    xlim([10,1015])
    grid on
    Nplots=15;
    heights=linspace(-2,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix0(heights*F.dh); % calculate round pixel planes
    for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    xticks(planes(1:2:end))
    xticklabels(ylab(1:2:end))
    xlabel('channel length [$d_h$]','Interpreter','latex')
    title('axial LFTP profiles')
    
    fname=sprintf('Fig%d axial profiles',fid);
    set(h,'PaperPosition', [1 1 12 9]);
    saveas(h,strcat(fname,'manual.png'))
    lft.f_BoFig2PDF(h,fname)
end

function lftmap(F)
    fid=307;h=figure(fid);clf;
    cas=2;
    botcut=8;
    upcut=1020;
    Nplots=15;
    heights=linspace(-2,12,Nplots); % at these channel heights we want to plot
    planes=lft.hcm2pix0(heights*F.dh); % calculate round pixel planes
    
    angb=linspace(-50,50,101); %angle basis
    angi=angb+90;             % angle index for array
    off=[0,0.5];
    for pin=[1,2]
        %ax=subplot('Position',[0.1+off(pin) 0.1 0.2 0.7]);
        ax=subplot(1,4,pin);
        im=squeeze(F.LFT(cas,botcut:upcut,pin,:));
        
        
        imshow(im,[]);
        set(gca,'YDir','normal');
        ylim([botcut,upcut]);
        caxis([0,120])

        if 1%pin==1
            xlabel('angular position [rad]')
        end
        xticks([1,91,181])
        xticklabels({'-\pi/2','0','\pi/2'})

        title(sprintf('Pin %d',pin))
        axis on
    
    
    
    for i=1:Nplots
        ylab{i}=num2str(heights(i));
    end
    
    if pin ==1
        yticks(planes(1:2:end))
        yticklabels(ylab(1:2:end))
        ylabel('channel length [$d_h$]','Interpreter','latex')
    else
        yticks([])
    end
    end
    
    cb=colorbar()
    ylabel(cb, 'LFT [$\mu m$]','Interpreter','latex')
    set(h,'Position',[0 0 320 660])
    fname=sprintf('Fig%d axial profiles',fid);
    set(h,'color','w');
    print(h,'filename','-dpdf','-r0')
    

%     u = h.Units;
%     h.Units = 'inches'
    export_fig Fig207_lftMap.png
    
end

function dh(F,d)
    % data:
    %d=f.loadSingleVariableMATFile('E:\20171213 final campaign\2_add_02.mat');
    im=squeeze(d(200:450,:,200));
    fid=201;
    h=figure(fid);clf;
    cas=2;
    height=[-1,0,4,8,12];
    planes=lft.hcm2pix0(F.dh*height);

    imshow(im',[])
    hold on
    set(gca,'YDir','normal')
    
    for i=1:length(planes)
        plot([0,size(d,1)],planes(i)*[1,1],'r')
    end
    axis on
    yticks(planes)
    yticklabels(height)
    ylabel('hydraulic diameters')
    xticks([])
    fname=sprintf('Fig%d pub horizontal profiles',fid);
    
    set(h,'color','w');
    saveas(h,strcat(fname,'manual4x4.png'))
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
        for cas=1:F.Ncases
           subplot(1,F.Ncases,cas)
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
           title(sprintf('case %d',cas))
           caxis(clim)
        end
        fname=sprintf('Fig%d LFT maps per pin %d',pin+26,pin);
        f.f_BoFig2PDF(h(cas),fname)
    end
    
    
end

function LFTOverview(F,T)
    fh=figure(73243);clf;
    height=20;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    
    xmarl=0.5;
    xmarr=0.2;
    xwid=1;
    xgap=0.1;
    xsum=sum([xmarl,xmarr,3*xwid,2*xgap]);
    spx=[xmarl,xmarl+xwid+xgap,xmarl+2*xwid+2*xgap]/xsum;
    
    ybot=0.2;
    yhig=1;
    ygap1=0.05;
    ycol=0.1;
    ygap2=0.2;
    ytop=0.1;
    ysum=sum([ybot,2*yhig+ygap1,ycol,ygap2,ytop]);
    spy=[ybot+yhig+ygap1+ycol+ygap2,ybot+ycol+ygap2]/ysum;
    
    dhs=[-3:3:9];
    planes=(lft.dh2plane(dhs)-10)/16;
    
    for cas=1:3
        for pin=1:2
            im=convn(squeeze(F.LFT(cas,:,pin,F.pth.angrang)),...
            mt.ElliKernelGen(1,5,1,2),'same');
            ax=subplot('position',[spx(cas),spy(pin),1/xsum,1/ysum]);
            imagesc(im,'AlphaData',~isnan(im));
            colormap jet
            caxis(max(F.LFT(:)*[0,1]))
            %colormap gray
            set(gca,'YDir','normal')
            hold on
            %plot([45 45],[1,63],'color',[1 1 1]*0.6)
            %plot([135 135],[1,63],'color',[1 1 1]*0.6)
            
            if cas==1
                yticks(planes)
                yticklabels(dhs)
                ylh=ylabel(sprintf('pin %d\n channel length [d_h]',pin));
                pub.BFylab(ylh,T.F)
            else
                yticks([]);
            end
            xticks([45 90 135])
            if pin==1
                %xticks([])
                
                xticklabels([])
                th=title(sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
                pub.BFtitle(th,T.F)
            else
                xlim();
                
                xticklabels({'-\pi/4','0','\pi/4'});
                xlh=xlabel('angular position');
                pub.BFxlab(xlh,T.F);
            end
            pub.BFaxis(ax,T.F);
            grid on
            ax.YGrid = 'off'
            
        end
    end
    %subplot('position',[xmarl/xsum,sum([ybot,yhig,ygap1])/ysum,1-(xmarl+xmarr)/xsum,ycol/ysum]);
    c=colorbar();
    c.Location='south';
    c.Position=[xmarl/xsum,ybot/ysum,1-(xmarl+xmarr)/xsum,ycol/ysum];
    ct=title(c,'Liqid film thickness [mm]');
    set(ct,'Units','data');
    pos = get (ct,'position');
    pos(2) = pos(2)-2.5;
    set(ct, 'position', pos);
    pub.BFxlab(ct,T.F)
    
    fname='LFTOverview';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    %open(sprintf('%s%s.pdf',T.F.saveto,fname))

end

function LFTlikeZboray(F,T)
    fh=figure(732431);clf;
    ax=gca();
    angs=(-45:45);
    
    dhs=[-3:3:9];
    planes=(lft.dh2plane(dhs)-10)/16;
    cas=1;
    pin=1;
    im1=squeeze(F.LFT(cas,:,pin,angs+90));
    pin=2;
    im2=squeeze(F.LFT(cas,:,pin,angs+90));
    
    im=convn([im2,im1],...
            mt.ElliKernelGen(1,5,1,2),'same');
    %im=im1;
    imagesc(im,'AlphaData',~isnan(im));
    colormap jet
    
    set(gca,'YDir','normal')
    hold on
    xticks([0 45 90 135 180])
    xticklabels(xticks+90)
    
    tx=[-10 0 20 40 60 80]
    yticks(lft.plane2rplane(lft.mm2plane(tx)));
    
    yticklabels({'-10','0','20','40','60','80'});
    
    xlabel('angle')
    ylabel('height [mm]')
    
    fname='LFTOverviewLikeZboray';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    %open(sprintf('%s%s.pdf',T.F.saveto,fname))

end

% pub.BFtitle(th,T.F)
% pub.BFfigure(fh,T.F)
% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)

function lftsurfplot(F,T,d,mtomo)
    % load in an empty tomo as demp(:,:,200:360)
    dvol=d(:,:,200:360);
    dz=16*0.127;
    r=F.geo.r_rod; % 5.14
    res1=0.127;         % resolution 1 or pixel width
    npix1=269;       % number of pixel withd of image
    numxmax1=(npix1-1)/2; % maximum pixel number in each direction
    numx1=linspace(-numxmax1,numxmax1,npix1); % number of each pixel
    coxmax1=res1*(npix1-1)/2; % coordinate of the furthest pixel
    cox1=linspace(-coxmax1,coxmax1,npix1); % coordinate list
    dvol = permute(dvol,[2 1 3]);
    for i=1:2;
        fh=figure(1296+i);clf
        height=20;
        width=8;
        set(fh,...
            'Units','centimeters',...
            'PaperPositionMode','manual',... %F.Fig.PPM,...
            'Position',[0 0 width height],...
            'PaperPosition',[0 0 width height],...
            'PaperSize',[width height]);
        
        [x,y,z]=meshgrid(cox1,cox1,((1:size(dvol,3))-92)*res1);
        p = patch(isosurface(x,y,z,dvol,0.095));
        isonormals(x,y,z,dvol,p);
        colormap jet
        p.FaceColor = 'red';
        p.EdgeColor = 'none';
        
    x2=x(150:end,150:end,:);
    y2=y(150:end,150:end,:);
    z2=z(150:end,150:end,:);
    dvol2=dvol(150:end,150:end,:);
    p2 = patch(isosurface(x2,y2,z2,dvol2,0.095));
    isonormals(x2,y2,z2,dvol2,p2);
    p2.FaceColor = 'green';
    p2.EdgeColor = 'none';
    daspect([1 1 1])
    view(3)
    axis tight
    camlight
    lighting gouraud
    xlh=xlabel('x [mm]');
    ylh=ylabel('y [mm]');
    if i==2
        ylh.Position=[-20 -3 -20];
    end
    zlh=zlabel('z [mm]');
    pub.BFylab(zlh,T.F);
    pub.BFxlab(xlh,T.F);
    pub.BFxlab(ylh,T.F);
    pub.BFaxis(gca,T.F);
    
    % rod surfaces
    cas=1;
    hold on
    cx=9.475*[1 0 -1 0];
    cy=9.475*[0 1 0 -1];
    falpha=[ 1 0.1; 0.1 1];
    for pin=i
        th=F.Pa.ProfAngles(pin,F.pth.angrang);
        z=((1:size(F.LFT,2))-F.par.spacerplane/16)*2;
        [Th Z R] = ndgrid(th, z,r);
        [XX,YY,ZZ] = pol2cart(Th,R,Z);
        XX=XX+cx(pin);
        YY=YY+cy(pin);
        surf(XX,YY,ZZ,permute(squeeze(F.LFT(cas,:,pin,F.pth.angrang)),[2,1]),...
            'edgecolor','none','facealpha',falpha(i,pin))
        caxis(max(F.LFT(~isnan(F.LFT(:)))*[0,1]))
    end
    %mesh(x, y, Z);
    xlim([-17,17])
    ylim(17*[-1 1])
    grid on
    than=title(sprintf('pin %d',i));
    pub.BFtitle(than,T.F)
    
    gx=cx(pin)+r*cos(F.Pa.ProfAngles(i,46));
    gy=cy(pin)+r*sin(F.Pa.ProfAngles(i,46));
    gz=[0,max(z)];
    plot3([gx gx],[gy gy],gz,'color',[0.6 0.6 0.6]);
    gx=cx(pin)+r*cos(F.Pa.ProfAngles(i,146));
    gy=cy(pin)+r*sin(F.Pa.ProfAngles(i,146));
  
    plot3([gx gx],[gy gy],gz,'color',[0.6 0.6 0.6]);
    if i==1
        view([-160,33.2])
        camlight
        view([-103,33.2])
    elseif i==2
        view([-165.5,33.2])
        camlight
        view([-11.5,33.2])
    end
    
    xc=linspace(-coxmax1,coxmax1,npix1);
    % add tomogram
    
     [xxx,yyy]=ndgrid(xc,xc);
     zzz=ones(npix1)*min(z(:));

     sp=surf(xxx,yyy,zzz,squeeze(d(:,:,300)),'edgecolor','none')
%     hold on
%     
%     % y-frame plane backside
%     yslice2=1;
%     surf(yx,yy+yslice2,yz,squeeze(d(:,yslice2,:)),'edgecolor','none')
    
% pub.BFfigure(fh,T.F)
% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)
    fname=sprintf('LFTOverview3d_%d',i);
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    end
    %open(sprintf('%s%s.pdf',T.F.saveto,fname))
end

function vanesplot(T,F,d,mtomo)
   
    dvol=d(:,:,250:360);
    dz=16*0.127;
    r=F.geo.r_rod; % 5.14
    res1=0.127;         % resolution 1 or pixel width
    npix1=269;       % number of pixel withd of image
    numxmax1=(npix1-1)/2; % maximum pixel number in each direction
    numx1=linspace(-numxmax1,numxmax1,npix1); % number of each pixel
    coxmax1=res1*(npix1-1)/2; % coordinate of the furthest pixel
    cox1=linspace(-coxmax1,coxmax1,npix1); % coordinate list
    dvol = permute(dvol,[2 1 3]);
    for i=1;
        fh=figure(11296+i);clf
%         height=20;
%         width=8;
%         set(fh,...
%             'Units','centimeters',...
%             'PaperPositionMode','manual',... %F.Fig.PPM,...
%             'Position',[0 0 width height],...
%             'PaperPosition',[0 0 width height],...
%             'PaperSize',[width height]);
        


%     
        [x,y,z]=meshgrid(cox1,cox1,((1:size(dvol,3)))*res1);
        p = patch(isosurface(x,y,z,dvol,0.095));
        isonormals(x,y,z,dvol,p);
        colormap jet
        p.FaceColor = 'red';
        p.EdgeColor = 'none';
        
    x2=x(150:end,150:end,:);
    y2=y(150:end,150:end,:);
    z2=z(150:end,150:end,:);
    dvol2=dvol(150:end,150:end,:);

    daspect([1 1 1])
    
    axis tight
    camlight
    lighting gouraud
%     xlh=xlabel('x [mm]');
%     ylh=ylabel('y [mm]');
    if i==2
        ylh.Position=[-20 -3 -20];
    end
    
    %xticks([])
    %yticks([])
    %zticks([])
    posx=7;
    posy=8;
    txh(1)=text(6.5,0,'pin 1','color',[1 1 1]);
    txh(2)=text(-2,posy,'pin 2','color',[1 1 1]);
    txh(3)=text(-10,0,'pin 3','color',[1 1 1]);
    txh(4)=text(-2,-posy,'pin 4','color',[1 1 1]);
                    
%     zlh=zlabel('z [mm]');
%     pub.BFylab(zlh,T.F);
%     pub.BFxlab(xlh,T.F);
%     pub.BFxlab(ylh,T.F);
    
    hold on
    lw=1;%linewidth
    l(1)=plot([-6,-4],[0 0],'color',[1 1 1]);%
    l(2)=plot([-7,-5],[2 4],'color',[1 1 1]);
    l(3)=plot([-7,-5],[-2 -4],'color',[1 1 1]);
    txh(5)=text(-3.5,0,'0°','color',[1 1 1]);
    txh(6)=text(-8,4.5,'\pi/4','color',[1 1 1]);
    txh(7)=text(-8,-4.5,'-\pi/4','color',[1 1 1]);
    pub.BFLine(l,T.F)
    
    pub.BFaxis(gca,T.F);
    set(txh,...
                    'FontUnits',T.Fig.yl.FU,...
                    'FontWeight',T.Fig.yl.FW,...
                    'FontSize',T.Fig.yl.FS-2,...
                    'FontName',T.Fig.yl.FN);
                 set(gcf, 'InvertHardCopy', 'off');
    cas=1;
    
   

    %mesh(x, y, Z);
    xlim([-17,17])
    ylim(17*[-1 1])
    grid on
    
    view([0,90])
    
    imshow((mtomo(:,:,2,25)-mtomo(:,:,1,25))',[],...
        'xdata',coxmax1*[-1,1],...
        'ydata',coxmax1*[-1,1])
    set(gca,'ydir','normal')
    
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
      
% pub.BFfigure(fh,T.F)
% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)
    fname=sprintf('/disc/vaneplot',i);
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    end
    
    
    
end


function PinCentered(F,T)
    fh=figure(12358);clf
    
    im=imread('V:\1phd\thesis\figures\disc\vaneplot_5.png');
    imshow(im)
    hold on
    lw=2;%linewidth
    
    l(1)=plot([165,185],[269/2 269/2],'color',[1 1 1]);%
     txh(1)=text(150,134,'0°','color',[1 1 1]);
     txh(5)=text(220,269/2,'pin 3','color',[1 1 1]);
     
    l(2)=plot([269/2,269/2],[165 185],'color',[1 1 1]);
     txh(2)=text(125,157,'270°','color',[1 1 1]);
     txh(6)=text(125,50,'pin 4','color',[1 1 1]);
     
    l(3)=plot([269/2,269/2],[85 105],'color',[1 1 1]);
     txh(3)=text(126,113,'90°','color',[1 1 1]);
     txh(7)=text(30,269/2,'pin 1','color',[1 1 1]);
     
    l(4)=plot([83,103],[269/2 269/2],'color',[1 1 1]);%
     txh(4)=text(105,134,'180°','color',[1 1 1]);
     txh(8)=text(125,220,'pin 2','color',[1 1 1]);

    pub.BFLine(l,T.F)
    
    pub.BFaxis(gca,T.F);
    set(txh,...
                    'FontUnits',T.Fig.yl.FU,...
                    'FontWeight',T.Fig.yl.FW,...
                    'FontSize',T.Fig.yl.FS-2,...
                    'FontName',T.Fig.yl.FN);
                 set(gcf, 'InvertHardCopy', 'off');
    
    
    
    
    
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
      
% pub.BFfigure(fh,T.F)
% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)
    fname=sprintf('/disc/vaneplotFinal',i);
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    
end

function ICONEtomos(F,T,mtomo,corr)
    % make some tomos for the presentation
    xr=31:239;
    
    im=[];
    plas=[7 22 32  58];
    for cas = 1:3
        imtemp=[];
        for pla = 1:length(plas)
            tom=(mtomo(xr,xr,cas+1,plas(pla))-mtomo(xr,xr,1,plas(pla)));
            imtemp=cat(1,imtemp,tom');
        end
        im=cat(2,im,imtemp);
    end
    fh=figure(1376);clf;
    imshow(-im,0.001*[-10 2.9]);set(gca,'YDir','normal')

    fpath='V:\1phd\conferences\ICONE26\tomos';
    savefig(fh,fpath)
    print(fpath,'-dpdf')
    print(fpath,'-dpng')
end
function ICONErad(F,T,mtomo,corr)
    
    % make figure with radio+ planes
    fh=figure(9655);
    plas=[7 22 32  58];
    imshow(exp(-corr(:,:,1))',[]);set(gca,'YDir','normal')
    hold on
    for pla=1:length(plas)
        plot([1 269], (plas(pla)*16-6)*[1 1],'r')
    end
    
    figure(2357);clf
    imagesc(squeeze(F.LFT(1,:,1,:)));set(gca,'YDir','normal')
    hold on
        for pla=1:length(plas)
        plot([1 269], plas(pla)*[1 1],'r')
        end
    xticks([])
    yticks([])
    
end

function lftsurfplot2(F,T,dvol)
    fh=figure(1297);clf
 
    dz=2;% roughly 16*0.127mm
    r=F.geo.r_rod; % 5.14
    res1=0.127;         % resolution 1 or pixel width
    npix1=269;
    pin=1;
    
    LFTdiff=F.LFT(2,:,:,:)-F.LFT(1,:,:,:);
    
    for pin=1:2
    th=F.Pa.ProfAngles(pin,:);
    z=(1:46)*2;
    
    
    [Th Z R] = ndgrid(th, z,r);
    
    [XX,YY,ZZ] = pol2cart(Th,R,Z);
    surf(XX,YY,ZZ,permute(squeeze(LFTdiff(1,:,1,:)),[2,1]),'edgecolor','none')
    caxis(mean(LFTdiff(:))*[0,5])
    hold on
    end
    %mesh(x, y, Z);
    xlim([-17,17])
    ylim(17*[-1 1])
    
end

function lftprofiles(F,T)
    fh=figure(5426); clf;
    height=20;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    
    plots=[5,22:6:size(F.pth.Profs,2)];
    ang1=-90:90;
    angu=ang1(F.pth.angrang);
    cols={'r','g','b'};
    
    xmarl=0.3;
    xwid=1;
    xgap=0.1;
    xmarr=0.2;
    xsum=xmarl+2*xwid+xgap+xmarr;
    
    ybot=1;
    yhig=1;
    ygap=0;
    ytop=1;
    ysum=ybot+length(plots)*yhig+(length(plots)-1)*ygap+ytop;
    
    xpos=[xmarl,xmarl+xwid+xgap]/xsum;
    ypos=zeros(1,length(plots));
    maxlft=0.16;%max(F.LFT(:));
    %maxlft=0.01;
    for i=1:length(plots)
        ypos(i)=(ybot+(i-1)*(yhig+ygap))/ysum;
    end
    for pin=1:2
        for plane=1:length(plots)
            ax=subplot('position',[xpos(pin),ypos(plane),xwid/xsum,yhig/ysum]);
            pub.BFaxis(ax,T.F)
            hold on
            for cas=1:3
                ph(cas)=plot(angu,squeeze(F.LFT(cas,plots(plane),pin,F.pth.angrang)),...
                    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
                pub.BFLine(ph(cas),T.F)
            end
            ylim([0,maxlft])
            xlim([min(angu),max(angu)])
            %plot([-45 -45],maxlft*[0 1],'color',0.6*[1 1 1])
            %plot([45 45],maxlft*[0 1],'color',0.6*[1 1 1])
            grid on
            xticks([-45 0 45])
            xticklabels([])
            if pin==1
                yticks([maxlft])
                if plane==floor(length(plots)/2)
                    ylh=ylabel('LFT [mm];    axial position ->');
                    pub.BFylab(ylh,T.F)
                end
                if plane==length(plots)
                    lh=legend(ph,'Location','best');
                    pub.BFlegend(lh,T.F)
                end
            else
                yticks([])
                txh=text(-70,maxlft*0.8,sprintf('%2.1f d_h',...
                    lft.plane2dh(plots(plane)*16+10)));
                set(txh,...
                    'FontUnits',T.Fig.yl.FU,...
                    'FontWeight',T.Fig.yl.FW,...
                    'FontSize',T.Fig.yl.FS-2,...
                    'FontName',T.Fig.yl.FN);
                
            end
            
            if plane==1
                xticks([-45 0 45]);
                xticklabels({'-\pi/4','0','\pi/4'});
                xlh=xlabel('angular position');
                pub.BFxlab(xlh,T.F);
            elseif plane==length(plots)
                th=title(sprintf('pin %d',pin));
                pub.BFtitle(th,T.F)
                
            end
        end
    end
    fname=sprintf('LFTangularprofiles',i);
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')

    
    
end

function lft360profiles(F,T)
    fh=figure(5489); clf;
    height=20;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    
    plots=[5,22:6:size(F.pth.Profs,2)];
    ang1=(-45:45)+90;
    angu=(-45:45)+90;
    angs=1:181;
    cols={'r','g','b'};
    
    xmarl=0.2;
    xwid=2;
    xgap=0.1;
    xmarr=0.1;
    xsum=xmarl+xwid+xmarr;
    
    ybot=0.5;
    yhig=1;
    ygap=0;
    ytop=0.5;
    ysum=ybot+length(plots)*yhig+(length(plots)-1)*ygap+ytop;
    
    xpos=[xmarl,xmarl+xwid+xgap]/xsum;
    ypos=zeros(1,length(plots));
    maxlft=0.16;%max(F.LFT(:));
    %maxlft=0.01;
    for i=1:length(plots)
        ypos(i)=(ybot+(i-1)*(yhig+ygap))/ysum;
    end
    
    for plane=1:length(plots)
        ax=subplot('position',[xpos(1),ypos(plane),xwid/xsum,yhig/ysum]);
        pub.BFaxis(ax,T.F)
        rectangle('Position',[45 0 90  maxlft],'EdgeColor','non',...
        'FaceColor',[0.9 0.9 1])
        rectangle('Position',[225 0 90  maxlft],'EdgeColor','non',...
        'FaceColor',[0.9 0.9 1])
        
        hold on
        
        for cas=1:3
            
            pin=3;
            ph(cas,1)=plot(0:45,squeeze(F.LFT(cas,plots(plane),pin,angs((45:90)+45))),...
                cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
            pin=4;
            ph(cas,2)=plot(45:135,squeeze(F.LFT(cas,plots(plane),pin,angs((0:90)+45))),...
                cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
            pin=1;
            ph(cas,3)=plot(135:225,squeeze(F.LFT(cas,plots(plane),pin,angs((0:90)+45))),...
                cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
            pin=2;
            ph(cas,4)=plot(225:315,squeeze(F.LFT(cas,plots(plane),pin,angs((0:90)+45))),...
                cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
            pin=3;
            ph(cas,5)=plot(315:360,squeeze(F.LFT(cas,plots(plane),pin,angs((0:45)+45))),...
                cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
            
            
            
            
        end
        pub.BFLine(ph(:),T.F)
        
        set(gca,'Layer','top')
        ylim([0,maxlft])
        xlim([0 360])
        
        grid on
        ax.XAxis.MinorTickValues = [45 135 225 315];
        ax.YAxis.MinorTickValues = [];
        grid minor
        
        %xticklabels([])
        yticks([maxlft])
        
        if plane==floor(length(plots)/2)
            ylh=ylabel('LFT [mm];    axial position ->');
            pub.BFylab(ylh,T.F)
        end
        if plane==length(plots)
            lh=legend(ph(:,1));
            pub.BFlegend(lh,T.F)
            set(lh,'Position',[0.5 0.8846 0.1703 0.0582])
            
        
            txp=-8    ;
            txh1(1)=text(15,maxlft*1.2,'pin 3');
            txh1(5)=text(330,maxlft*1.2,'pin 3');
            txh1(2)=text(txp+90,maxlft*1.2,'pin 4');
            txh1(3)=text(txp+180,maxlft*1.2,'pin 1');
            txh1(4)=text(txp+270,maxlft*1.2,'pin 2');
            set(txh1,...
                'FontUnits',T.Fig.ti.FU,...
                'FontWeight','bold',...
                'FontSize',T.Fig.ti.FS,...
                'FontName',T.Fig.ti.FN);
        end
        txh=text(10,...
            1.2*max(squeeze(F.LFT(1,plots(plane),pin,angs((45:90)+45)))),...
            sprintf('%2.1f d_h',...
            lft.plane2dh(plots(plane)*16+10)));
        set(txh,...
            'FontUnits',T.Fig.yl.FU,...
            'FontWeight',T.Fig.yl.FW,...
            'FontSize',T.Fig.yl.FS-2,...
            'FontName',T.Fig.yl.FN);
        
        xticks([0 90 180 270 360])
        if plane==1
            
           ;
            xlh=xlabel('angular position [°]');
            pub.BFxlab(xlh,T.F);
        else
           xticklabels([]);
            
        end
    end
    fname=sprintf('disc/LFT360profiles',i);
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')

    
    
end

function ZDD(T,F,C)
% combined Zboray Damsohn Dotox profile plot
%read in Zboray stuff
z=xlsread('G:\cbolesch\Data\Zboray_cold_neutrons_LFT_profile.xlsx');
% experiment 3
 zdatax=[z(:,13);nan;z(:,15)];
 zdatay=[z(:,14);nan;z(:,16)];

dhd=18.85; % hydraulic diameter double subchannel
dh=F.geo.dh*10;

angs=1:181;
ang=-45:45;
offset=0;

% looking for zboray's area from 60-100mm in DH's
zheights=[60 100]/dhd; % 3.1830    5.3050
[zrang]=round(lft.plane2rplane(lft.dh2plane(zheights)));
zrange=zrang(1):zrang(2);



fh=figure(7894);clf;
pub.BFfigure(fh,T.F)
height=6;
width=16;
 set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);

ax=subplot(1,2,1);
hold on

th=title(sprintf('%1.1f-%1.1f Dh downstream',zheights(1),zheights(2)));
pub.BFtitle(th,T.F)

cols={'r','g','b'};
for cas=1:3
p(cas)=plot(ang+45+90+45,squeeze(mean(mean(F.LFT(cas,zrange,[1,3],ang+90),3),2)),...
    cols{cas},'displayname',sprintf('case',cas));
q(cas)=plot(ang+45+45,squeeze(mean(mean(F.LFT(cas,30:40,[2,4],ang+90),3),2)),...
    cols{cas});
end
r=plot(zdatax-45,zdatay/1000,'displayname','cold neutrons');

% gas=1;
% jd=5; % j_g=20 m/s
% bet=5; % beta=0.004
% r(2)=plot(C.angles-45,flip(C.exp(6,gas,jd,bet).mean(:,20)/1000,1));
% r(3)=plot(C.angles-45,flip(C.exp(6,gas,jd+1,bet+1).mean(:,20)/1000,1));


% % insert calvin stuff
% x=C.angles(9:16)+45;
% y1=mean(C.exp(1,5,5).mean(9:16,(1:10)+10)/1000,2);
% y2=mean(C.exp(1,5,9).mean(9:16,(1:10)+10)/1000,2);
% y=0.75*y1+0.25*y2;
% % plot(x,y,...
% %     '--','color','k','displayname','calvin');


xlh=xlabel('angle [deg]');
xticks([45 90 135 180 225 ])
xlim([45 225]);
%lh=legend([p,r]);
grid on   

ylh=ylabel('LFT [mm]');
ylim([0 0.4])

%pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)
pub.BFLine([p,q,r],T.F)


% no spacer
ax=subplot(1,2,2);

th=title('upstream');
pub.BFtitle(th,T.F)
zdx=[z(:,17);nan;z(:,19)];
zdy=[z(:,18);nan;z(:,20)];
hold on
ylim([0 0.4])
cols={'r','g','b'};
zrange=2:9;
angs=1:181;
ppin=[-1 0 1 2];
for cas=1:3
   
    pin=3;
    p(cas,1)=plot(0:45,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:45)+45)),2)),...
    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=4;
    p(cas,2)=plot(45:135,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2)),...
    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=1;
    p(cas,3)=plot(135:225,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2)),...
    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=2;
    p(cas,4)=plot(225:315,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2)),...
    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=3;
    p(cas,5)=plot(315:360,squeeze(mean(F.LFT(cas,zrange,pin,angs((45:90)+45)),2)),...
    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    



end
% p(cas)=plot(ang+45+90+90,squeeze(mean(mean(F.LFT(cas,zrange,[1,3],ang+90),3),2)),...
%     cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
% q(cas)=plot(ang+45+90,squeeze(mean(mean(F.LFT(cas,zrange,[2,4],ang+90),3),2)),...
%     cols{cas});

r(1)=plot(zdx-45,zdy/1000,'displayname','Zboray','color',[0    0.4470    0.7410]);
r(2)=plot(zdx(zdx<45)-45+360,zdy(zdx<45)/1000,'displayname','Zboray','color',[0    0.4470    0.7410]);
% r(3)=plot(C.angles-45,flip(C.exp(1,gas,jd,bet).mean(:,20)/1000,1));
% r(4)=plot(C.angles-45,flip(C.exp(1,gas,jd+1,bet+1).mean(:,20)/1000,1));
xlim([0 360])
xlh=xlabel('angle [deg]');
xticks([0 90 180 270 360])
%xlim([90 270]);
lh=legend([p(:,1)',r(1)]);
grid on   

ylh=ylabel('LFT [mm]');
 pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)
pub.BFLine([p(:)',q,r],T.F)
 

for i=1:2
    
    ax=subplot(1,2,i);
    pos=ax.Position;
    pos(2)=0.2;
    pos(4)=0.7;
    
    set(ax,'Position',pos)
    
end


 fname='lit/LFTZobray';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


    
end

function Calvin(F,T,C)
    angs=1:181;
    ang=-45:45;
    offset=0;
    
    dhd=18.85; % hydraulic diameter double subchannel
    dh=F.geo.dh*10;

    % looking for zboray's area from 60-100mm in DH's
    zheights=[60 100]/dhd; % 3.1830    5.3050
    [zrang]=round(lft.plane2rplane(lft.dh2plane(zheights)));
    zrange=zrang(1):zrang(2);
    zundis=2:8;
    
    fh=figure(78941);clf;
    pub.BFfigure(fh,T.F)
    height=3*6;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    
    ymax=0.03;
    cols={'r','g','b'};
    liq={'air','C4F8','He'};
    cstr={'Re_{gl}: Calvin','We_l: Calvin','We_l: Calvin'};
    Dbetas=[0.0069,0.0047,0.0030];
    
    zCstart=round(find(C.heights>60,1,'first')/120*64);
    zCstop=round(find(C.heights<100,1,'last')/120*64);
    xrang1=9:16;
    xrang2=1:16;
    xrang=1:8;
    
    spacs=[6,6,6];
    js=[4 7 7]; % 20    35
    bs=[9 5 3]; % 20    40    80
    gas=[1 1 3];
    subs=[1 3 5];
    for i=1:3 % Re  or We scling
        ax=subplot(3,2,subs(i));
        hold on
        grid on
        
        for cas=1:3
            p(cas)=plot(ang+45+90+45-180,squeeze(mean(mean(F.LFT(cas,zrange,[1,3],ang+90),3),2))/dh,...
                cols{cas},'displayname',sprintf('case',cas));
            q(cas)=plot(ang+45+45,squeeze(mean(mean(F.LFT(cas,30:40,[2,4],ang+90),3),2))/dh,...
                cols{cas});
        end
        
        for b=1:3
            r(b)=plot(C.angles(xrang1)-45,...
                flip(mean(...
                C.exp(spacs(i),gas(i),js(i),bs(b)).mean(xrang,[zCstart:zCstop])/1000/C.Dh,2)),...
                '--','color',cols{b});
            s(b)=plot(C.angles(1:8)-45,...
                mean(...
                C.exp(spacs(i),gas(i),js(i),bs(b)).mean(9:16,[zCstart:zCstop])/1000/C.Dh,2),...
                '--','color',cols{b},'displayname',...
                sprintf('Calvin %0.3f',C.exp(spacs(i),gas(i),js(i),bs(b)).beta/10000));
        end
        %r(1)=plot(zdatax-45,zdatay/1000/dhz,'displayname','cold neutrons');
        
        
        
        if i==1
            th=title('with spacer');
            ylh=ylabel({sprintf('%s %s J_g= %2.0f m/s',...
                cstr{i},liq{gas(i)},C.exp(spacs(i),gas(i),js(i),bs(i)).jg),...
                'LFT [D_h]'});
            pub.BFylab(ylh,T.F)
        end
        
        if any(i==[2,3]);
            
            ylh=ylabel({sprintf('%s %s J_g= %2.0f m/s',...
                cstr{i},liq{gas(i)},C.exp(spacs(i),gas(i),js(i),bs(i)).jg),...
                'LFT [D_h]'});
            pub.BFylab(ylh,T.F)
        end
        
        if i==3
            xlh=xlabel('angle [deg]');
            pub.BFxlab(xlh,T.F)
        end
        %
        xticks([-45 0 45 90 135 180 225 ])
        xlim([-45 135]);
        ylim([0 ymax])
        
        if i==2;
            lh(i)=legend([s]);pub.BFlegend(lh(i),T.F)
        end
        %set(lh(i),'Position',[0.3185 0.85 0.11 0.06])
        
        pub.BFtitle(th,T.F)
        pub.BFaxis(ax,T.F)
        pub.BFLine([p,q,r,s],T.F)
        
        ax=subplot(3,2,subs(i)+1);
        hold on
        grid on
        
        if i==1
            th=title('without spacer');
        end
        
        
        for cas=1:3
            p(cas)=plot(ang+45+90+45-180,squeeze(mean(mean(F.LFT(cas,zundis,[1,3],ang+90),3),2))/dh,...
                cols{cas},'displayname',sprintf('DoToX %.4f',Dbetas(cas)));
            q(cas)=plot(ang+45+45,squeeze(mean(mean(F.LFT(cas,zundis,[2,4],ang+90),3),2))/dh,...
                cols{cas});
        end
        
        spa=1;
        for b=1:3
            r(b)=plot(C.angles(xrang2)-45,...
                flip(mean(...
                C.exp(1,gas(i),js(i),bs(b)).mean(xrang2,[zCstart:zCstop])/1000/C.Dh,2)),...
                '--','color',cols{b},'displayname',...
                sprintf('%d %d %d',gas(i),js(i),bs(b)));

        end
        
        xticks([-45 0 45 90 135 180 225 ])
        xlim([-45 135]);
        if i==3
            xlh=xlabel('angle [deg]');
            pub.BFxlab(xlh,T.F)
        end
        
        grid on
        ylim([0 ymax])

        %lh=legend([r]); pub.BFlegend(lh,T.F)
        pub.BFaxis(ax,T.F)
        
        pub.BFLine([p,q,r],T.F)
        if i==2;
            lh(i)=legend([p]);pub.BFlegend(lh(i),T.F)
        end
        
    end
    
%     for i=2
%         set(lh(i),'Position',get(lh(i),'Position')+[0.0 0.02 0 0])
%     end
    
    fname='lit/Calvin';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end
%
function ZDDscaled(T,F,C)
    % combined Zboray Damsohn Dotox profile plot, scaled
    %read in Zboray stuff
    z=xlsread('G:\cbolesch\Data\Zboray_cold_neutrons_LFT_profile.xlsx');
    % experiment 3
    zdatax=[z(:,13);nan;z(:,15)];
    zdatay=[z(:,14);nan;z(:,16)];
    
    dhratio=0.44936; % infinite lattice
    %dhratio=0.0.4841; % actuual channel
    dhz=18.85;
    dh=8.47;
    
    angs=1:181;
    ang=-45:45;
    offset=0;
    % looking for zboray's area from 60-100mm in DH's
    dhd=18.85;
    zheights=[60 100]/dhd; % 3.1830    5.3050
    [zrang]=round(lft.plane2rplane(lft.dh2plane(zheights)));
    zrange=zrang(1):zrang(2);
    
    fh=figure(78953);clf;
    pub.BFfigure(fh,T.F)
    height=6;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    
    ax=subplot(1,2,1);
    hold on
    
    th=title(sprintf('%1.1f-%1.1f Dh downstream',zheights(1),zheights(2)));
    pub.BFtitle(th,T.F)
    ymax=0.02;
    
    cols={'r','g','b'};
    for cas=1:3
        p(cas)=plot(ang+45+90+45,squeeze(mean(mean(F.LFT(cas,zrange,[1,3],ang+90),3),2))/dh,...
            cols{cas},'displayname',sprintf('case',cas));
        q(cas)=plot(ang+45+45,squeeze(mean(mean(F.LFT(cas,30:40,[2,4],ang+90),3),2))/dh,...
            cols{cas});
end
r(1)=plot(zdatax-45,zdatay/1000/dhz,'displayname','cold neutrons');
gas=1;
jd=8; 
bet=5;
spa=6
h=20;
mm=h/64*120;

xlh=xlabel('angle [deg]');
xticks([45 90 135 180 225 ])
xlim([45 225]);
%lh=legend([p,r]);
grid on  
ylim([0 ymax])

ylh=ylabel('LFT [D_h]');

%pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)
pub.BFLine([p,q,r],T.F)

y1=flip(C.exp(spa,gas,jd,bet).mean(:,h)/1000/dhz,1);
y2=flip(C.exp(spa,gas,jd,9).mean(:,h)/1000/dhz,1);
w1=1;
w2=0;
y=w1*y1+w2*y2;
% 
%  r(2)=plot(C.angles-45,y,'displayname',sprintf('%d m/s, %0.4f @%2.1fD_h',C.exp(spa,gas,jd,bet).jg,...
%         C.exp(spa,gas,jd,bet).beta/10000,mm/18.85));
%  pub.BFLine([p,q,r],T.F)
%  legend(r(2))
%  y1=flip(C.exp(6,gas,jd,bet).mean(:,20)/1000/dhz,1);
% y2=flip(C.exp(6,gas,jd,9).mean(:,20)/1000/dhz,1);
% y=0.75*y1+0.25*y2;
%  r(3)=plot(C.angles-45,y);
%  pub.BFLine([p,q,r],T.F)

% no spacer
ax=subplot(1,2,2);

th=title('upstream');
pub.BFtitle(th,T.F)
zdx=[z(:,17);nan;z(:,19)];
zdy=[z(:,18);nan;z(:,20)];
hold on
ylim([0 ymax])
cols={'r','g','b'};
zrange=2:9;
angs=1:181;
ppin=[-1 0 1 2];
for cas=1:3
    
    pin=3;
    p(cas,1)=plot(0:45,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:45)+45)),2))/dh,...
        cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=4;
    p(cas,2)=plot(45:135,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2))/dh,...
        cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=1;
    p(cas,3)=plot(135:225,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2))/dh,...
        cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=2;
    p(cas,4)=plot(225:315,squeeze(mean(F.LFT(cas,zrange,pin,angs((0:90)+45)),2))/dh,...
        cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    pin=3;
    p(cas,5)=plot(315:360,squeeze(mean(F.LFT(cas,zrange,pin,angs((45:90)+45)),2))/dh,...
        cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
    
end
% p(cas)=plot(ang+45+90+90,squeeze(mean(mean(F.LFT(cas,zrange,[1,3],ang+90),3),2)),...
%     cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
% q(cas)=plot(ang+45+90,squeeze(mean(mean(F.LFT(cas,zrange,[2,4],ang+90),3),2)),...
%     cols{cas});

r(1)=plot(zdx-45,zdy/1000/dhz,'displayname','Zboray','color',[0    0.4470    0.7410]);
r(2)=plot(zdx(zdx<45)-45+360,zdy(zdx<45)/1000/dhz,'displayname','Zboray','color',[0    0.4470    0.7410]);

% y1=flip(C.exp(1,gas,jd,bet).mean(:,40)/1000/dhz,1);
% y2=flip(C.exp(1,gas,jd,9).mean(:,40)/1000/dhz,1);
% y=w1*y1+w2*y2;
% 
%  r(3)=plot(C.angles-45,y,'displayname',sprintf('%d m/s, %0.4f',C.exp(spa,gas,jd,bet).jg,...
%         C.exp(spa,gas,jd,bet).beta/10000));

%r(3)=plot(C.angles-45,flip(C.exp(1,gas,jd,bet).mean(:,40)/1000/dhz,1));

xlim([0 360])
xlh=xlabel('angle [deg]');
xticks([0 90 180 270 360])
%xlim([90 270]);
lh=legend([p(:,1)',r([1])],'location','northeast');
set(lh,'Position',[0.7135 0.62 0.1702 0.2533])
grid on

ylh=ylabel('LFT [D_h]');
pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)
pub.BFLine([p(:)',q,r],T.F)


for i=1:2
    
    ax=subplot(1,2,i);
    pos=ax.Position;
    pos(2)=0.2;
    pos(4)=0.7;
    
    set(ax,'Position',pos)
    
end


fname='lit/LFTZobrayScaled';
fpath=sprintf('%s%s',T.F.saveto,fname);

savefig(fh,fpath)
print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
print(sprintf('%s%s',T.F.saveto,fname),'-dpng')



end

function CalvinNoSpacer(F,T,C,Z)
fh=figure(2781);clf
ax=gca;
hold on
grid on
js=2:12;
cols=hsv(length(js));
liqs={'air','C4F8','He'};
b=6;
gas=1;
spa=1;

for j=1:length(js)
    q(j)=plot(C.angles,mean(C.exp(spa,gas,js(j),b).mean,2)','o-','color',cols(j,:),...
        'displayname',sprintf('%d m/s, %0.4f',C.exp(spa,gas,js(j),b).jg,...
        C.exp(spa,gas,js(j),b).beta/10000));
end
q(j+1)=plot(Z.ux,Z.uy,'k','displayname','zboray 26.6m/s 0.0050');

legend(q)
title(sprintf('scan of j_g without spacer, %s',liqs{gas}))
xlabel('angle')
ylabel('mean LFT')

 fname='cal/nospacer_scan_j';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


fh=figure(2782);clf
ax=gca;
hold on
grid on
js=2:12;
bets=2:14;
cols=hsv(length(bets));
j=5;

for b=1:length(bets)
    p(b)=plot(C.angles,mean(C.exp(spa,gas,j,bets(b)).mean,2)','o-','color',cols(b,:),...
        'displayname',sprintf('%d m/s, %0.4f',C.exp(spa,gas,j,bets(b)).jg,...
        C.exp(spa,gas,j,bets(b)).beta/10000));
end
p(b+1)=plot(Z.ux,Z.uy,'k','displayname','zboray 26.6m/s 0.0050');
legend(p)
title(sprintf('scan of 1-\beta without spacer, %s',liqs{gas}))
xlabel('angle')
ylabel('mean LFT')

 fname='cal/nospacer_scan_b';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


figure(2783);clf
ax=gca;
hold on
grid on
heights=[5 10 20 30 40 50 60 ];
mm=heights/64*120;

j=8;
b=5;
spa=6;

nh=length(heights);
cols=hsv(nh);
for h=1:nh
    r(h)=plot(C.angles,flip(C.exp(spa,gas,j,b).mean(:,heights(h))),'color',cols(h,:),...
        'displayname',sprintf('%d m/s, %0.4f @%2.0fmm',C.exp(spa,gas,j,b).jg,...
        C.exp(spa,gas,j,b).beta/10000,mm(h)));
end

r(h+1)=plot(Z.x,Z.y,'k','displayname','zboray 26.6m/s 0.0050');
legend(r)
title(sprintf('%d m/s, %0.4f,%s with spacer',C.exp(1,1,j,b).jg,...
    C.exp(1,1,j,b).beta/10000,liqs{gas}))
xlabel('angle')
ylabel('mean LFT')

 fname='cal/spacerZvC';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    
    
    figure(2784);clf
ax=gca;
hold on
grid on
heights=[5 10 20 30 40 50 60 ];
mm=heights/64*120;

j1=6;
b1=6;
w1=2/3;

j2=3;
b2=6;
w2=1/3;

spa=6;
nh=length(heights);
cols=hsv(nh);
r=[];
for h=1:nh
    y1=flip(C.exp(spa,gas,j1,b1).mean(:,heights(h)));
    y2=flip(C.exp(spa,gas,j2,b2).mean(:,heights(h)));
    y=w1*y1+w2*y2;
    jg=w1*C.exp(spa,gas,j1,b1).jg+w2*C.exp(spa,gas,j2,b2).jg;
    beta=(w1*C.exp(spa,gas,j1,b1).beta+w2*C.exp(spa,gas,j2,b2).beta)/10000;
    r(h)=plot(C.angles,y,'color',cols(h,:),...
        'displayname',sprintf('%d m/s, %0.4f @%2.0fmm',jg,beta,mm(h)));
end

r(h+1)=plot(Z.x,Z.y,'k','displayname','zboray 26.6m/s 0.0050');
legend(r)
title({sprintf('%d m/s, %0.4f, %s with spacer',jg,beta,liqs{gas}),...
    'interpolated from 15 and 30 m/s',})
xlabel('angle')
ylabel('mean LFT')

 fname='cal/spacerZvC2';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


end

function CloneFig(inFigNum,OutFigNum)
% this program copies a figure to another figure
% example: CloneFig(1,4) would copy Fig. 1 to Fig. 4
% Matt Fetterman, 2009
% pretty much taken from Matlab Technical solutions:
% http://www.mathworks.com/support/solutions/en/data/1-1UTBOL/?solution=1-1UTBOL
hf1=figure(inFigNum);
hf2=figure(OutFigNum);
clf;
pub.compCopy(hf1,hf2);
end

function compCopy(op, np)
%COMPCOPY copies a figure object represented by "op" and its % descendants to another figure "np" preserving the same hierarchy.

ch = get(op, 'children');
if ~isempty(ch)
nh = copyobj(ch,np);
for k = 1:length(ch)
pub.compCopy(ch(k),nh(k));
end
end;
return;
end



function lftprofiles2(F,T)
    fh=figure(5427); clf;
    height=20;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);

     plots=[5,22:6:size(F.pth.Profs,2)];
    ang1=-90:90;
    angu=ang1(F.pth.angrang);
    cols={'r','g','b'};
    
    xmarl=0.3;
    xwid=1;
    xgap=0.1;
    xmarr=0.2;
    xsum=xmarl+2*xwid+xgap+xmarr;
    
    ybot=1;
    yhig=1;
    ygap=0;
    ytop=1;
    ysum=ybot+length(plots)*yhig+(length(plots)-1)*ygap+ytop;
    
    xpos=[xmarl,xmarl+xwid+xgap]/xsum;
    ypos=zeros(1,length(plots));
    maxlft=max(F.LFT(:));
    for i=1:length(plots)
        ypos(i)=(ybot+(i-1)*(yhig+ygap))/ysum;
    end
    
    
    for pin=1:2
        
        for plane=1:length(plots)
            ax=subplot('position',[xpos(pin),ypos(plane),xwid/xsum,yhig/ysum]);
            pub.BFaxis(ax,T.F)
            hold on
            for cas=1:3
                ph(cas)=plot(angu,squeeze(F.LFT(cas,plots(plane),pin,F.pth.angrang)),...
                    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas)));
                ph(cas)=plot(angu,squeeze(F.LFT(cas,plots(plane),pin+2,F.pth.angrang)),...
                    cols{cas},'displayname',sprintf('%3.2f kg/min',F.par.flowkgs(cas)));
            end
            ylim([0,maxlft])
            xlim([min(angu),max(angu)])
            %plot([-45 -45],maxlft*[0 1],'color',0.6*[1 1 1])
            %plot([45 45],maxlft*[0 1],'color',0.6*[1 1 1])
            grid on
            xticks([-45 0 45])
            xticklabels([])
            if pin==1
                yticks([0,maxlft])
                if plane==floor(length(plots)/2)
                    ylh=ylabel('LFT [mm];    axial position ->');
                    pub.BFylab(ylh,T.F)
                end
                if plane==length(plots)
                    lh=legend(ph,'Location','best');
                    pub.BFlegend(lh,T.F)
                end
            else
                yticks([])
                text(-70,maxlft*0.8,sprintf('%d mm',(plots(plane)-6)*2))
                
            end
            
            if plane==1
                xticks([-45 0 45]);
                xticklabels({'-\pi/4','0','\pi/4'});
                xlh=xlabel('angular position');
                pub.BFxlab(xlh,T.F);
            elseif plane==length(plots)
                th=title(sprintf('pin %d',pin));
                pub.BFtitle(th,T.F)
                
            end
        end
    end
    title('the profiles are symmetric')
    % pub.BFtitle(th,T.F)
    % pub.BFfigure(fh,T.F)
    % pub.BFlegend(lh,T.F)
    % pub.BFaxis(ax,T.F)
    % pub.BFylab(ylh,T.F)
    % pub.BFxlab(xlh,T.F)
    % pub.BFLine(p,T.F)
    
    
end

function [fh,ax]=checkprofile(im,profs,x,y)
    figure(1347);clf;
    s=size(im,1);
    [xx,yy]=ndgrid(1:s,1:s);
    zz=zeros(s,s);
    su=surf(xx,yy,zz,im,'edgecolor','none');
    colormap(gray)

    nang=size(profs,3);
    col=hsv(nang*4);
    hold on
    for pin=1:4
        for ang2=1:nang;
            plot3(squeeze(x(:,pin,ang2)),squeeze(y(:,pin,ang2)),...
                squeeze(profs(:,pin,ang2)),...
                'color',col(ang2+(pin-1)*nang,:))
        end
    end
    fh=gcf;
    ax=gca;
end


function volumePlot(d,T)
    % before this function execute the following:
    % load('T')
%     T.cas=6;
%     T.rep=2;
%     d=f.ReadXrayTomo(T);
%     perm=[3,1,2];
%     dperm=permute(d(:,31:1015,85:(85+1357)),perm);
    
    fid=50;
    fig=figure(fid);clf;
    pub.BFfigure(fig,T.F)
    %set(gcf,'pos',[80   103   600   668])
    
%     sorig=size(d);
    perm=[3,1,2];
%     dperm=permute(d,perm);
    s=size(d);
    planes=[300,50,1357];
    planes=[planes(perm(1)),planes(perm(2)),planes(perm(3))];
    
    colormap gray
    s=size(d);
    
    % y-frame plane 
    [yx,yz]=ndgrid(1:s(1),1:s(3));
    yy=zeros(s(1),s(3));
    yslice=280;
    surf(yx,yy+yslice,yz,squeeze(d(:,yslice,:)),'edgecolor','none')
    hold on
    
    % y-frame plane backside
    yslice2=1;
    surf(yx,yy+yslice2,yz,squeeze(d(:,yslice2,:)),'edgecolor','none')

    % back radio
    xplane=planes(1);
    [xy,xz]=ndgrid(1:s(2),1:s(3));
    xx=zeros(s(2),s(3));
    surf(xx+xplane,xy,xz,squeeze(d(xplane,:,:)),'edgecolor','none')
    
    % front radio
    range=1:460;
    surf(xx(range,:),xy(range,:),xz(range,:)+1,...
        squeeze(d(1,range,:)),'edgecolor','none') 

    % front radio, lower slice
    yrange=1:100;
    surf(xx(:,yrange),xy(:,yrange),xz(:,yrange)+1,...
        squeeze(d(1,:,yrange)),'edgecolor','none') 
    
    % bottom sinogram
    zplane=1;
    [zx,zy]=ndgrid(1:s(1),1:s(2));
    zz=zeros(s(1),s(2));
    surf(zx,zy,zz+zplane,squeeze(d(:,:,zplane)),'edgecolor','none')
    
        
    % middle sinogram
    xrange2=190:460   ;
    for zplane3=[270,600,985];
    surf(zx(:,xrange2),zy(:,xrange2),zz(:,xrange2)+zplane3,...
        squeeze(d(:,xrange2,zplane3)),'edgecolor','none')
    end
    
   % [frames, x, y]
   % total frame
   pub.cube_plot(gca,[0 0 0 ],[1357,640,985],[0,0,0],2)
   % spacer
   pub.cube_plot(gca,[0 200 130 ],[1357,250,180],[0,0,1],1)
   % angular gauge
   pub.cube_plot(gca,[0 0 0 ],[1357,640,60],[0,1,0],1)
   % gauge
   pub.cube_plot(gca,[0 0 300 ],[1357,180,460],[1,0,1],1)
   
    xlim([0,1357])
    ylim([0,640])
    zlim([0,985])

    xticks([])
    yticks([])
    zticks([])

    set(gca,'YDir','reverse')
    set(gca,'FontSize',14)
    axis equal
    view(-58.3000,21.2000)
    
    tightInset = get(gca, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(gca, 'Position', position);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    

   

   print(sprintf('%sdatacube.pdf',T.Fig.saveto),'-dpdf')
   open(sprintf('%sdatacube.pdf',T.Fig.saveto))
end

function VectorGraphicTest
    fig=figure(1);clf;
    set(fig,...
        'renderer','opengl')
    k=1000;
    data =rand(k);
    range=1:k;
    [xx,yy]=ndgrid(range,range);
    zz=zeros(size(xx));1
    surf(xx,yy,zz,data,'EdgeColor','none')
    hold on
    surf(xx,yy,zz+1000,data,'EdgeColor','none')
    surf(xx,yy,zz+500,data,'EdgeColor','none')
    print('TestVecIMage.pdf','-dpdf')
    print('TestVecIMage.png','-dpng')
    open('TestVecIMage.pdf')
end

function PubPlot(h,ax,fpath)
    set(h,'FontSize',14)
    tightInset = get(ax, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(ax, 'Position', position);
    set(ax,'units','centimeters')
    pos = get(ax,'Position');
    ti = get(ax,'TightInset');
    
    set(h, 'PaperUnits','centimeters');
    set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
   h.Renderer='Painters';
   
   print(h,fpath,'-dpdf')
   print(h,fpath,'-dpng')
   savefig(h,fpath)
end

function spcerVid(F,T,corr)
    % make a video of the spacer
    corr=corr(:,6:end-5,:);
    fh=figure(26235)
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 269 size(corr,2)];
    fh.Units='pixel';
    fh.Position=[50 50 269 size(corr,2)];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    crange=[0 1]
    len=size(corr,3);
    M(len) = struct('cdata',[],'colormap',[]);
    
    v = VideoWriter('V:\1phd\conferences\ICONE26\spacervid');
    v.FrameRate=30;
        open(v)
    for j = 1:50%len
        imshow(corr(:,:,j)',[crange])
        set(gca,'YDir','normal')
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame)
    end
    close(v)
end

function cube_plot(ax,orig,dims,col,lw)
% plots a cube starting at orig with extent dims
% cube_plot(ax,orig,dims,color,lw)
% ax = axes handle
% orig = [x,y,z] start point
% dims = [x,y,z] sizes
% col = [r g b] color
% lw = line width

X=orig(1);
Y=orig(2);
Z=orig(3);

x=dims(1);
y=dims(2);
z=dims(3);

% lower square
plot3([X,X+x,X+x,X,X],[Y,Y,Y+y,Y+y,Y],[Z,Z,Z,Z,Z],'color',col,'linewidth',lw)
% upper square
plot3([X,X+x,X+x,X,X],[Y,Y,Y+y,Y+y,Y],[Z+z,Z+z,Z+z,Z+z,Z+z],'color',col,'linewidth',lw)
% legs
plot3([X,X],[Y,Y],[Z,Z+z],'color',col,'linewidth',lw)
plot3([X+x,X+x],[Y,Y],[Z,Z+z],'color',col,'linewidth',lw)
plot3([X,X],[Y+y,Y+y],[Z,Z+z],'color',col,'linewidth',lw)
plot3([X+x,X+x],[Y+y,Y+y],[Z,Z+z],'color',col,'linewidth',lw)

end

function printVector
    fig=figure();
    
    % plot an image, ready for 3d plot
    img = imread('ngc6543a.jpg');
    s=size(img);
    [xx,yy]=ndgrid(1:s(1),1:s(2));
    zz=zeros(s(1:2));
    %surf(xx,yy,zz,img,'edgecolor','none')
    
    
    % add 3d content
    x=rand(10,1)*s(1);
    y=rand(10,1)*s(2);
    z=rand(10,1)*s(1);
    plot3(x,y,z)
    hold on
    imshow(img)
    
    %export as vector graphic
    fig1.Renderer='Painters';
    print(fig,'vectorgraphic','-dpdf')
end
    
function BeamNDoseNormalisation(T)
    cas=1;
    rep=6;
    
    spd=0.02; % sub plot delta

    h=openfig('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\Fig 007 01 06 dose area checks.fig');
    axesObjs = get(h, 'Children');
    dataObjs = get(axesObjs, 'Children');
    fid=1000;
    fh=figure(fid);clf;
    pub.BFfigure(fh,T.F)
    
        
    % dose vs frame plot
    s1=subplot('Position',[0.11,0.1,0.8,0.2]);
    p=plot(dataObjs{2}.XData,dataObjs{2}.YData/mean(dataObjs{2}.YData));
    pub.BFLine(p,T.F)
    
    grid on
    xlim([1,T.Raw.nframes(cas,rep)]);
    xticks([1 500 1000]);
    title('Dose variation during one measurement')
    xlh=xlabel('frames');
    pub.BFxlab(xlh,T.F)
    
    ylh=ylabel('rel. dose');
    pub.BFylab(ylh,T.F)
    close(h);
    pub.BFaxis(gca,T.F)

    % surf plot
    h=openfig('C:\Users\cbolesch\Desktop\HS14_calibration desktop\final campaign\Fig 006 01 06 Beam fit.fig');
    axesObjs = get(h, 'Children');
    dataObjs = get(axesObjs, 'Children');
       
    figure(fid);
    ax2=subplot('Position',[0.1,0.4,0.45,0.5]);
    
    mask=T.dose.mask; %re-label
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];
    rf =200; % reduction factor
    im=single(T.Raw.rad(:,:,cas,rep));
    im(mask==0)=[];
    [x3,y3]=ndgrid([1,20:20:640],[1,32:32:1024]);
    dose=feval(T.fit.fit{cas,rep}{1},x3,y3);
    meandose=mean(dose(:));
    
    s1=surf(x3,y3,dose/meandose,'DisplayName','fit');
    hold on
    sc1=plot3(xx(1:rf:end)',yy(1:rf:end)',im(1:rf:end)'/meandose,'ob','DisplayName',...
        'fit base points');
    %p1=plot(T.fit.fit{1,6}{1},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)');
    grid on
    view([22,48])
    lh=legend([s1,sc1]);
    set(lh,'Location','northeast');
    pub.BFlegend(lh,T.F)
    
        close(h)
    %surf properties
    set(s1,'EdgeColor','none',...
              'FaceAlpha',0.6)
    set(sc1,'MarkerFaceColor',[0 0 1],...
            'MarkerEdgecolor',[1 1 1],...
            'MarkerSize',6)      
    % axes properties    
    set(ax2,'xticklabel',{[]},...
        'yticklabel',{[]},...
        'TickDir','in')
    pub.BFaxis(gca,T.F)
    
    th=title('Beam non-uniformity fit');
    pub.BFtitle(th,T.F)
    
    a=single(T.Raw.rad(:,:,cas,rep));
    ylim([1,1024])
    xlim([1,620])
   
    xlh=xlabel('x');
    set(xlh,'Position',[400,-50 0.9]);
    pub.BFxlab(xlh,T.F)
    
    ylh=ylabel('y');
    set(ylh,'Position',[670,450,0.9]);
    pub.BFxlab(xlh,T.F)
    
    zlh=zlabel('local relative beam intensity');
    pub.BFylab(zlh,T.F)
    
    % radiography with rectangles
    s3=subplot('Position',[0.605,0.4,0.35,0.5]);
    p2=imshow(squeeze(T.Raw.rad(:,:,cas,rep))',[]);
    set(s3,'xticklabel',{[]},...
        'yticklabel',{[]},...
        'YDir','normal');
    th=title('Open Beam Areas');
    pub.BFtitle(th,T.F)
    
    for i=1:3
        r=rectangle('Position',[T.dose.windows(1,i),...
            T.dose.windows(3,i),...
            T.dose.windows(2,i)-T.dose.windows(1,i),...
            T.dose.windows(4,i)-T.dose.windows(3,i)],...
            'Facecolor','b');
    end

    xlh=xlabel('x');
    set(xlh,'Position',[320.5003 0 1]);
    pub.BFxlab(xlh,T.F)
    
    ylh=ylabel('y');
    pub.BFylab(ylh,T.F)

    pub.BFaxis(gca,T.F)
    
    print(sprintf('%sDoseNBeam.eps',T.Fig.saveto),'-depsc2')
    %print(sprintf('%sDoseNBeam.pdf',T.Fig.saveto),'-dpdf')
    print(sprintf('%sDoseNBeam.pdf',T.Fig.saveto),'-dpdf')
    %open(sprintf('%sDoseNBeam.pdf',T.Fig.saveto))

    
end

function shiftplot(T)
    % 
    fid=1001
    h=figure(fid),clf
    pub.BFfigure(h,T.F)
    set(h,'Position',[0 0 2*T.F.Fig.FigW T.F.Fig.FigH],...
        'PaperPosition',[0 0 2*T.F.Fig.FigW T.F.Fig.FigH])
    cas=6;T.cas=cas;
    rep=2;T.rep=rep;

    
    sino=squeeze(T.Raw.rad(:,T.Raw.CropRange(T.cas,T.rep,1,2):...
        T.Raw.CropRange(T.cas,T.rep,2,2),cas,rep));
    nframes=size(sino,2);
    
    sinor=squeeze(T.Raw.rad(:,T.Raw.CropRange(1,1,1,2):...
        T.Raw.CropRange(1,1,2,2),1,1));
    nframesr=sizer(sino,2);
    
    subplot(1,3,1:2)
    imshow(sino(1:end,80:1437),[])
    
    xlh=xlabel('angle');
    pub.BFxlab(xlh,T.F)
    
    ylh=xlabel('x position');
    pub.BFylab(ylh,T.F)
    
    title('two sinograms with offset')
    
    
end

function QualityMask(rec,dif,mask,T)
   % makes a plot of the quality mask
   tstring={'parallel','Fan'};
   for i=1:size(mask,3);
       disp(i)
       fid=1001;
       fh=figure(fid);clf;
       pub.BFfigure(fh,T.F);
       fheight=9;
       fh.Position(4)=fheight;
       fh.PaperPosition(4)=fheight;
       fh.PaperSize(2)=fheight;
       
       subplot('position',[0.025,0.1,0.3,0.8  ]);
       imshow(rec(:,:,i)',[]);set(gca,'YDir','normal');
       
       xlh=xlabel('reconstruction');
       pub.BFxlab(xlh,T.F);
       %set(xlh,'Position',[400,-50 0.9]);
       
       subplot('position',[0.35,0.1,0.3,0.8 ]);
       imshow(dif(:,:,i)',[]);set(gca,'YDir','normal');
       xlh=xlabel('gradient');
       pub.BFxlab(xlh,T.F);
       title(sprintf('%s',tstring{i}))
       
       subplot('position',[0.675,0.1,0.3,0.8  ]);
       imshow(mask(:,:,i)',[]);set(gca,'YDir','normal');
       xlh=xlabel('edge quality mask');
       pub.BFxlab(xlh,T.F);
       
       print(sprintf('%sQualitymask%s.pdf',T.F.saveto,tstring{i}),'-dpdf')
       open(sprintf('%sQualitymask%s.pdf',T.F.saveto,tstring{i}))
   end
   
end

function [ax,fh]=decenter(T,d)
    % parforhacks
    sinsh=T.Rec.SinoRotShift;
    ang=deg2rad(T.Rec.angles);
    recsize=T.Rec.recsize;
    detPitch=T.Rec.detPitch;
    src=T.Rec.src;
    det=T.Rec.det;
    niter=T.Cen.niter;
    BS=T.Raw.BS;
    xmin=T.Cen.xmin;
    xmax=T.Cen.xmax;
    MaskThresh=T.Cen.MaskThresh;
    fh=figure(236);clf
    
    
    shift=[0,0.1,3];
    im=[];
    
    fh.Units='pixel';
    fh.PaperUnits='centimeters';
    papsize=2.8;
    fh.Position=[0 0 recsize recsize];
    fh.PaperPosition=[0 0 papsize papsize];
    fh.PaperSize=[papsize papsize];
    
    for i=1:length(shift)
        rec(:,:,i)=a.FBPexplFan(...
            f.fraccircshift(squeeze(d(:,400,:)),shift(i))',...
            recsize,ang,detPitch,src,det);
        imshow(rec(:,:,i),[-1 10]*1e-06)
        set(gca,'YDir','normal')
        ax=gca;
        ax.Units='pixel';
        ax.Position=[0 0 recsize recsize];
        fh.Units='pixel';
        fh.Position=[50 50 recsize recsize];
        fh.Units='centimeters';
        fh.PaperPosition=fh.Position;
        fh.PaperSize=fh.Position(3:4);
%         fh.InnerPosition=fh.Position;
%         fh.OuterPosition=fh.Position;
%          set(fh, 'PaperUnits', 'normalized')
%          set(fh, 'PaperPosition', [0 0 1 1])
         set(fh, 'PaperUnits', 'centimeters')
         set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])

        fname=sprintf('shift%d',i);
        print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
        open(sprintf('%s%s.pdf',T.F.saveto,fname))
        
    end
end

function voidprofile(F,T,mtomo)
    fh=figure(4283);clf;
    height=9;
    width=16;
    set(fh,...
        'Units','centimeters',...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 width height],...
        'PaperPosition',[0 0 width height],...
        'PaperSize',[width height]);
    col=hsv(4);
    ph=[];
    dhs=-3:3:9;
    planes=(lft.dh2plane(dhs)-10)/16;
    for cas=1:3
        ax=subplot(1,3,cas);
        pub.BFaxis(ax,T.F);
        hold on
        ymax=max(F.void.mean(:))/F.XS*100;
        for reg=1:4
           
            ph(reg)=plot(squeeze(F.void.mean(reg,cas,:)+F.void.mean(reg+4,cas,:))/2/F.XS*100,'color',col(reg,:),...
                'Displayname',sprintf('reg. %d',reg));
            
        end
        %sum
        ph(5)=plot(squeeze(mean(F.void.mean(:,cas,:),1))/F.XS*100,':','color',[0 0 0],...
                'Displayname','channel mean');
            pub.BFLine(ph,T.F)
        ph(6)=area([10 21], ymax*[1 1],'Facecolor',0.9*[1 1 1],...
            'edgecolor','none','displayname','spacer');    
        
%         rectangle('Position',[10,0.41,12,max(F.void.mean(:))/0.058*100],...
%             'Facecolor',0.9*[1 1 1],'Edgecolor','none')
        ylim([0,ymax]);
        xlim([1,63]);
        xticks(planes);
        xticklabels(dhs);
        grid on
        th=title(sprintf('%3.2f kg/min',F.par.flowkgs(cas+1)));
        pub.BFtitle(th,T.F)
        if cas==1
        ylh=ylabel('liquid fraction [%]');
        pub.BFylab(ylh,T.F)
        elseif cas==2
            xlh=xlabel('axial position [d_h]');
            pub.BFxlab(xlh,T.F)
        end
        if cas==1
            lh=legend(ph,'location','southwest');
        end
        set(gca,'Layer','top')
    end
    
    subplot(1,3,1)
    %lh=legend(ph,'location','southwest');
    pub.BFlegend(lh,T.F);
    fname='VoidProfile';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    
    
    % make a geometrical map for the regions
    regiomap=repmat(mtomo(:,:,2,34)-mtomo(:,:,1,34),[1,1,3]);
    regiomap=(regiomap-min(regiomap(:)))/(max(regiomap(:))-min(regiomap(:)));
    figure(8328);clf
    imshow(permute(regiomap,[2 1 3]))
    set(gca,'YDir','normal')
    
    colo=hsv(4);
    for c=1:3
        for reg=1:4
            im=zeros(269);
            im=or((F.void.finalmask==reg),(F.void.finalmask==reg+4))*colo(reg,c);
            
            regiomap(:,:,c)=regiomap(:,:,c)+im;
        end
    end
    regiomap=(regiomap-min(regiomap(:)))/(max(regiomap(:))-min(regiomap(:)));
    
    fh=figure(8329);clf
    %     height=9;
    %     width=16;
    %     set(fh,...
    %         'Units','centimeters',...
    %         'PaperPositionMode','manual',... %F.Fig.PPM,...
    %         'Position',[0 0 width height],...
    %         'PaperPosition',[0 0 width height],...
    %         'PaperSize',[width height]);
    imshow(permute(regiomap,[2 1 3]),[])
    set(gca,'YDir','normal')
    text(185,151, '1','Fontsize',14,'color',colo(1,:))
    text(150,184, '2','Fontsize',14,'color',colo(2,:))
    text(120,184, '3','Fontsize',14,'color',colo(3,:))
    text(80,151, '4','Fontsize',14,'color',colo(4,:))
    
    
    
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    
    fname='VoidProfileMap';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end

function misalignment(F,T,im)
        fh=figure(8330);clf
    %     height=9;
    %     width=16;
    %     set(fh,...
    %         'Units','centimeters',...
    %         'PaperPositionMode','manual',... %F.Fig.PPM,...
    %         'Position',[0 0 width height],...
    %         'PaperPosition',[0 0 width height],...
    %         'PaperSize',[width height]);
    imshow(im',[])
    set(gca,'YDir','normal')
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    
    fname='misalignment';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end

function radiograph(T,F,im)
        fh=figure(8331);clf
    %     height=9;
    %     width=16;
    %     set(fh,...
    %         'Units','centimeters',...
    %         'PaperPositionMode','manual',... %F.Fig.PPM,...
    %         'Position',[0 0 width height],...
    %         'PaperPosition',[0 0 width height],...
    %         'PaperSize',[width height]);
    imshow(im',[])
    set(gca,'YDir','normal')
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    
    
    fname='radiograph';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end


function paths1(F,T,mtomo)
    % show one path
    [xx,yy]=ndgrid(1:F.recsize,1:F.recsize);
    zz=zeros(F.recsize,F.recsize);
    fh=figure(1365);clf;

cas=4;
pla=22;
mtpla=pla;
dcas=1;

%su=surf(xx,yy,zz,mtomo(:,:,cas,pla),'edgecolor','none');
ax=imshow( mtomo(:,:,2,mtpla)'-mtomo(:,:,1,mtpla)',[]);
set(gca,'YDir','normal')
ax2=gca();
clim2=ax2.CLim;

%colormap(gray)
angs=11:10:F.pth.n_angles-10;
nang=length(angs);
profang=10;
col=hsv(nang);

hold on
for pin=3
for ang2=1:nang;
    ang=angs(ang2);
    x=linspace(F.cen.cen(pin,1),F.Pa.ProfEndPnt(pin,ang,1),size(F.pth.Profs,5));
    y=linspace(F.cen.cen(pin,2),F.Pa.ProfEndPnt(pin,ang,2),size(F.pth.Profs,5));
 %   plot3(x,y,squeeze(F.difProf(dcas,pla,pin,ang,:)),'color',col(ang2,:))
    p3=plot(x,y,'color',[0 0 1]);
    pub.BFLine(p3,T.F)
    if ang2==profang
        refang=ang; % for later
        p4=plot(x,y,'color',[1 0 0]);
        pub.BFLine(p4,T.F)
    end
end
  %squeeze(F.difProf(dcas,pla,pin,ang,:))
 %squeeze(F.difProf(dcas,pla,pin,ang,:))
end
%ylim([0 150])
    recsize=F.recsize;
    ax=gca;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    
    fname='Paths1';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


if 0    
% show the profiles
fh=figure(1366);clf
 pub.BFfigure2(fh,T.F,1)
      fh.Position(1)=10;

ax=subplot(1,2,1);
pub.BFaxis(ax,T.F)
hold on
p(1)=plot((1:size(F.difProf,5))*F.recres,squeeze(F.pth.Profs(2,pla,pin,refang,:)),...
    'red','DisplayName','w/  film');
p(2)=plot((1:size(F.difProf,5))*F.recres,squeeze(F.pth.Profs(1,pla,pin,refang,:)),...
    ':','DisplayName','w/o film');


pub.BFLine(p,T.F)
ylh=ylabel('\mu [1/mm]');
pub.BFylab(ylh,T.F)
xlh=xlabel('path length [mm]');
pub.BFxlab(xlh,T.F)

grid on
ylm=ax.YLim;
ylim(ylm+[-0.001 0])
lh=legend('location','best');
pub.BFlegend(lh,T.F)
legend('boxoff')

ax=subplot(1,2,2);
pub.BFaxis(ax,T.F)
p2=plot((1:size(F.difProf,5))*F.recres,squeeze(F.difProf(1,pla,pin,refang,:)),...
    'DisplayName','film only');
pub.BFLine(p2,T.F)
lh=xlabel('path length [mm]');
pub.BFxlab(xlh,T.F)
grid on
ylim(ylm+[-0.001 ]);
yticklabels([]);

lh=legend('location','best');
pub.BFlegend(lh,T.F)
legend('boxoff')

    fname='Paths2';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end
    

% show the profiles
fh=figure(1367);clf
 pub.BFfigure2(fh,T.F,1.15)
      fh.Position(1)=10;

ax=gca()
pub.BFaxis(ax,T.F)
hold on
p(1)=plot((1:size(F.difProf,5))*F.recres,squeeze(F.difProf(cas-1,pla,pin,refang,:)),...
    'red','DisplayName','path profile');
areax=F.pth.lftSumRange;
areay=max(F.pth.lftSumRange):(F.pth.r_path+1);
area((areax)*F.recres,squeeze(F.difProf(cas-1,pla,pin,refang,areax)),...
    'DisplayName','film','facealpha',0.5,'edgecolor','none');
area((areay)*F.recres,squeeze(F.difProf(cas-1,pla,pin,refang,areay)),...
    'DisplayName','gas core','facealpha',0.5,'edgecolor','none');

pub.BFLine(p,T.F)
ylh=ylabel('\mu [1/mm]');
pub.BFylab(ylh,T.F)
xlh=xlabel('distance from pin center [mm]');
pub.BFxlab(xlh,T.F)

grid on
ylm=ax.YLim;
%ylim(ylm+[-0.001 0])
lh=legend('location','northwest');
yylim=ax.YLim;

patch([4.14 5.14 5.14 4.14],[yylim(1) yylim(1) yylim(2) yylim(2)],...
    0.5*[1 1 1],'edgecolor','none', 'facealpha',0.5,'displayname','wall');
pub.BFlegend(lh,T.F)
    fname='Paths3';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end

function sumrange(F,T,mtomo) % scan for good LFT summing range
    fh=figure(1367);clf
    pub.BFfigure2(fh,T.F,3);
    fh.Position=[-50.8000 0 50.8000 26.5642];
    %fh.Position(1)=10;
    anglist=1:20:181;
    anglist([1,end])=[];

    nang=length(anglist);
    areax=F.pth.lftSumRange;  % film area
    areay=(max(areax)+1):(F.pth.r_path+1);% void area
    pla=8;
    pin=1;
    for cas=1:3
        for ang =1:nang
            nplot=(cas-1)*nang+ang;
            ax=subplot(3,nang,nplot);
            plot((1:size(F.difProf,5))*F.recres,...
                squeeze(F.difProf(cas,pla,pin,anglist(ang),:)),...
                'red','DisplayName','path profile');
            hold on
            area((areax)*F.recres,squeeze(F.difProf(cas,pla,pin,anglist(ang),areax)),...
                'DisplayName','film','facealpha',0.5,'edgecolor','none');
            area((areay)*F.recres,squeeze(F.difProf(cas,pla,pin,anglist(ang),areay)),...
                'DisplayName','void');
            ylim(1e-3*[-2 4])
            xlim([2 8])
            grid on
            title(sprintf('cas %d ang %d°',cas,anglist(ang)))
            
            ax.XTick=sort([ax.XTick(1),ax.XTick(end),5.14,4.14]);
            
            
        end
    end
    
    ax=gca();
    pub.BFaxis(ax,T.F)
    hold on
    
    ylh=ylabel('\mu [1/mm]');
    pub.BFylab(ylh,T.F)
    xlh=xlabel('path length [mm]');
    pub.BFxlab(xlh,T.F)
    
    grid on
    ylm=ax.YLim;
    %ylim(ylm+[-0.001 0])
    lh=legend('location','northwest');
    % add tomogram
    figure(1368);clf
    for cas=1:3
        subplot(1,3,cas)
    imshow(mtomo(:,:,cas+1,pla)'-mtomo(:,:,1,pla)',[])
    set(gca,'YDir','normal')
    title(sprintf('plane %d, pin %d',pla,pin))
    end
    
    
end

function ZborayVoid(T,F) 
fh=figure(2647);clf
ax=gca;
hold on

col=[1 0 0;
      0 1 0;
      0 0 1];

for cas=1:3
    p(cas)=plot(F.ivoid.ivoid(cas,:),lft.plane2mm(lft.rplane2plane(1:F.h)),...
        '--','color',col(cas,:),...
        'displayname',sprintf('%3.2f kg/min e',F.par.flowkgs(cas+1)));
    p(cas+3)=plot(squeeze(mean(F.void.mean(:,cas,:),1))/F.XS,...
    lft.plane2mm(lft.rplane2plane(1:F.h)),...
    'color',col(cas,:),'displayname',sprintf('%3.2f kg/min e_g',F.par.flowkgs(cas+1)));
end
grid on
ylim([-30 90])
lh=legend('location','southeast');
set(lh,'Color',0.9*[1 1 1])

xlim([0 0.12])
xticks([])
yticks([])
% yticks(ytx);
% yticklabels(ytx)
 pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
 pub.BFLine(p,T.F)
set(gcf, 'InvertHardCopy', 'off')

fname='lit/ZborayVoidss';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')



end

function im=ImSegExp(im,n)
%expands the image n imes in all directions
for i=1:n
im2=circshift(im,-1,1);
im3=circshift(im,1,1);
im4=circshift(im,-1,2);
im5=circshift(im,1,2);
im=(im+im2+im3+im4+im5)>0;
end
im=double(im);

end

function ivoid
end




function RobertTomo(F,T,mtomo)
    % take the figures from robert and turn them intopublication grade
    % plots
    
    %1: new six-tomo

     
     figfilename='G:\cbolesch\Data\201712xx_final_campaign\r2magic_factor\fig2.fig';
     fig1=openfig(figfilename,'invisible');
         fh=figure(2364);
     pub.BFfigure2(fh,T.F,2)
     imtot=[];
     ex=25; %extracrop
     for i=1:5
         im=imrotate(fig1.Children(i).Children.CData,45,'bilinear','crop');
         resize=269;
         im=imresize(im,[resize,resize]);
         diff=(F.recsize-resize)/2;
         if mod(diff,2)==1
             error('make resize odd')
         end
         im=im(1+diff:end-diff,1+diff:end-diff);
         imtot(:,:,6-i)=im(1+ex:end-ex,1+ex:end-ex);
     end
     close(fig1)
     %add our own tomo
     pla=9;
     cas=2;
     imtot(:,:,6)=mtomo(1+ex:end-ex,1+ex:end-ex,cas,pla)-...
         mtomo(1+ex:end-ex,1+ex:end-ex,1,pla);
     
     %merge images
     immerge=[];
     
     for i=1:2
         immergetemp=[];
        for j=1:3
            k=j+(i-1)*3;
            
           % size(immergetemp)
            immergetemp=cat(2,immergetemp,imtot(:,:,k));
            %size(immergetemp)
        end
       % size(immerge)
        immerge=cat(1,immerge,immergetemp);
        %size(immerge)
     end
    
     fh=figure(2364);
         pub.BFfigure2(fh,T.F,2)
     
     imshow(immerge,[]);
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[0 0 16 16*2/3];
    fh.Position=[0 0 16 16*2/3];
    fh.PaperPosition=[0 0 16 16*2/3];

    
    
     hold on
     xx=[30 400];
     lfts={'0.1 mm','0.2 mm','0.3 mm','0.4 mm','0.5 mm','measurement'};
    for i = 1:3
        for j=1:2
            k=i+(j-1)*3;
            if k==6
                continue
            end
    t=text(size(imtot,1)*(i-1+0.5)-25,xx(j),lfts{k},'color',[1 1 1])   ;
    pub.BFxlab(t,T.F)
        end
    end
    set(gcf, 'InvertHardCopy', 'off');
    
    fname='validationtomo';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)
    
end
function robertMu(T,F,mtomo)
    % plot spectrum, detector response, and the mu
    figfilename='G:\cbolesch\Data\201712xx_final_campaign\r2magic_factor\fig1.fig';
    fig1=openfig(figfilename,'invisible');
    
    lines=fig1.Children(8).findobj; %0.2mm plot
    specx=lines(2).XData;
    specy=lines(2).YData;
    effx=lines(3).XData;
    effy=lines(3).YData;
    mu3=lines(4).YData(1);
    mu2=lines(5).YData(1);
    mu1=lines(6).YData(1);
    xsx=lines(7).XData;
    xsy=lines(7).YData;
    close(fig1)
    
    %evil plot hack
    %specy(specy==max(specy))=2;
    
    fh=figure(2364);clf;
    pub.BFfigure2(fh,T.F,1)
    ax=gca()
    grid on
    hold on
    xlim([20 110])
    ylim([0 3])
    
    cols=get(groot,'DefaultAxesColorOrder');
    p(1)=plot(specx,specy/1000000/2.5,'color',cols(3,:),'displayname', 'relatve x-ray spectrum');
    p(2)=plot(effx,effy,'color',cols(1,:),'displayname','detector efficiency');
    p(3)=plot(xsx,xsy,'color',cols(2,:),'displayname','chloroform attenuation');
    p(4)=plot([20 110],mu2*[1 1],':','color',cols(2,:),'displayname','literature value')
    p(5)=plot([20 110],mu3*[1 1],'--','color',cols(2,:),'displayname','simulation result')
    pub.BFLine(p,T.F)
    
    annotation('arrow',0.95*[1 1],[0.5 0.41]);
    annotation('arrow',0.95*[1 1],[0.24 0.33]);

    t=text(71,1.3,{'setup specific','beam hardening'});
    pub.BFxlab(t,T.F)
    
    xlh=xlabel('energy [keV]')
    ylh=ylabel('attenuation [1/cm]; efficiency [-]');
    yticks([0 1 2 3])
    lh=legend(p,'location','northeast');
    pub.BFlegend(lh,T.F)
    pub.BFaxis(ax,T.F)
    pub.BFylab(ylh,T.F)
    pub.BFxlab(xlh,T.F)
    ax.Position=[ax.Position(1:2),0.98-ax.Position(1:2)];
    
    fname='validationsimulation';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
    
end

function robertMu2(T,F,mtomo)
    % the mu versus lft plot
    figfilename='G:\cbolesch\Data\201712xx_final_campaign\r2magic_factor\fig6.fig';
    fig1=openfig(figfilename,'invisible');
    
    mu1x=fig1.Children(2).Children(1).XData;
    mu1y=fig1.Children(2).Children(1).YData;
    
    mu2x=fig1.Children(2).Children(2).XData;
    mu2y=fig1.Children(2).Children(2).YData;
    
    fh=figure(2365);clf;
    pub.BFfigure2(fh,T.F,1)
    ax=gca()
    grid on
    hold on
    
    p(2)=plot(mu1x,mu1y,'x-','Displayname','simulation')
    hold on
    p(1)=plot(mu2x,mu2y,'Displayname','literature')
    ylim([0 1]);
    pub.BFLine(p,T.F)
     xlh=xlabel('film thickness [mm]');
    ylh=ylabel('attenuation [1/cm]');

    lh=legend(p,'location','southeast');
    pub.BFlegend(lh,T.F)
    pub.BFaxis(ax,T.F)
    pub.BFylab(ylh,T.F)
    pub.BFxlab(xlh,T.F)
    ax.Position=[ax.Position(1:2),0.97-ax.Position(1:2)];
    
        fname='validationMu';
    fpath=sprintf('%s%s',T.F.saveto,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s',T.F.saveto,fname),'-dpng')
end

function CADAlignment(F,T,mtomo,P2)
    fh=figure(1368);clf;

cas=2;
pla=6;

%su=surf(xx,yy,zz,mtomo(:,:,cas,pla),'edgecolor','none');
imshow( mtomo(:,:,cas,pla)',[])
set(gca,'YDir','normal')
hold on
    
%plot centers
for i =F.cen.pcenid%1:length(P2.c)
    plot(P2.c(i,1),F.cen.P2.c(i,2),'ob');
    %text(F.cen.P2.c(i,1),F.cen.P2.c(i,2),num2str(F.cen.P2.r(i)),'FontSize',20);
end

f42_PlotPhan(P2,gca,[1 0 0],[1 4]);

end


function DataHist(d,T,F)
% make a histogram of the data
fh=figure(5648);clf;
ax=gca;
edges=linspace(0,2^16,101);
histogram(d(200:500,400:1000,:),edges,'Normalization','probability',...
    'displayname','channel region')
hold on
histogram(d(200:500,100:1000,:),edges,'Normalization','probability',...
    'displayname','channel + spacer')
lh=legend('location','northwest');
xlh=xlabel('counts [#]');
ylh=ylabel('probability');
th=title('distribution of detector counts distribution');
grid on
xlim([0, 2^16])
pub.BFtitle(th,T.F)
%pub.BFfigure2(fh,T.F,2)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)
pub.BFlegend(lh,T.F)
x=2;
y=1;
    set(fh,...
        'Units',T.F.Fig.Units,...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 x*T.F.Fig.FigW y*T.F.Fig.FigH],...
        'PaperPosition',[0 0 x*T.F.Fig.FigW y*T.F.Fig.FigH],...
        'PaperSize',[x*8 y*6]);
    fname='DataHist';
    fpath=sprintf('%s%s%s',T.F.saveto,'technical\',fname);
    % the figure is 1.9 GB large - forbidden to save .fig
    % savefig(fh,fpath)
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpdf')
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpng')
end

function im=CA(im) % Contrast Adjust
    %im=(im-min(im(:)))/(max(im(:))-min(im(:)));
    
end

function NeutronXray(T,F,mtomo)
    
    fh=figure(7861);clf
    pub.BFfigure2(fh,T.F,2)
    rs=F.recsize;
    %empty=imresize(imread('V:\1phd\thesis\figures\neutrons\empty_vs_opb_SLICE_0600.png'),[rs rs]);
    %cl3d2o=imresize(imread('V:\1phd\thesis\figures\neutrons\CHCL3_vs_d2o_SLICE_0600.png'),[rs rs]);
    %cl3ob=imresize(imread('V:\1phd\thesis\figures\neutrons\CHCL3_vs_opb_SLICE_0600.png'),[rs rs]);
    
    % on new setup
    empty=imread('G:\cbolesch\HS14_calibration desktop\final campaign\NED paper\empty_vs_opb_SLICE_0600.png');
    cl3d2o=imread('G:\cbolesch\HS14_calibration desktop\final campaign\NED paper\CHCL3_vs_d2o_SLICE_0600.png');
    cl3ob=imread('G:\cbolesch\HS14_calibration desktop\final campaign\NED paper\CHCL3_vs_opb_SLICE_0600.png');
    
    zoom=0.36;
    empty=imresize(empty,zoom);
    cl3d2o=imresize(cl3d2o,zoom);
    cl3ob=imresize(cl3ob,zoom);
    
    angle=58;
    empty=imrotate(empty,angle,'bilinear');
    cl3d2o=imrotate(cl3d2o,angle,'bilinear');
    cl3ob=imrotate(cl3ob,angle,'bilinear');
    
    imsize=size(empty);
    
    cropbox=[floor((imsize(1)-rs)/2),floor((imsize(2)-rs)/2),rs-1,rs-1];
    
    shift=6;
    empty=circshift(imcrop(empty,cropbox),shift,2);
    cl3d2o=circshift(imcrop(cl3d2o,cropbox),shift,2);
    cl3ob=circshift(imcrop(cl3ob,cropbox),shift,2);
    
    %disp(size(empty))
    
    T.F.saveto = 'G:\cbolesch\HS14_calibration desktop\final campaign\NED paper\';
    
    im1=mtomo(:,:,2,8);
    im2=mtomo(:,:,1,8);
    im3=im1-im2;
    
    im=cat(1,cat(2,empty,cl3ob,cl3d2o),100000*9*cat(2,im1,im2,5*im3));
    ax=gca;
    imshow(im,[]);
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[1 1 16-2 16*2/3-2];
    fh.Position=[0 0 16 16*2/3];
    fh.PaperPosition=[0 0 16 16*2/3];
    xticks(rs*[0.5 1.5 2.5])
    xticklabels({'w/o chloroform','with chloroform','difference'});
    t=annotation('textarrow',[0.26 0.29],[0.87 0.8],...
        'String','displacement plug','color',[1 1 1]);
    yticks(rs*[0.5 1.5]);
    yticklabels({'neutrons','xrays'})
    ytickangle(90)
    axis on
    set(gca,'TickLength',[0 0])
    


pub.BFaxis(ax,T.F)

    
    
        fname='NeutVsXray';
        folder='';
    fpath=sprintf('%s%s%s',T.F.saveto,folder,fname);
    savefig(fh,fpath)
    print(sprintf('%s%s%s',T.F.saveto,folder,fname),'-dpdf')
    print(sprintf('%s%s%s',T.F.saveto,folder,fname),'-dpng')
end

function dynamicbias(T)
    % creates the dynamic bias plot
    
    fh=figure(38467);clf
    pub.BFfigure2(fh,T.F,1)
    ax=gca;
    pub.BFaxis(ax,T.F)
    hold on
    
    n=20;
    
    of=0.5;
    mu=0.8;
    d=1;
    dens=n/d;
    X=0.8;
    ymax=1.2;
    xof=linspace(d,d+of,round(dens*of));
    yof=ones(size(xof));
    
    
        rectangle('Position',[0 0 X*d  ymax],'EdgeColor','non',...
        'FaceColor',[0.8 0.8 1])
    x=linspace(0,X*d,n);
    x2=linspace(0,d,n);
    y2=exp(-mu*x2);
    
    y3=X*y2+(1-X);
    
    y=exp(-mu*x);
    p(1)=plot(x,y,'b','Displayname','instant image');
    p(5)=plot([-of,0],[1,1],'b');
    p(6)=plot([X*d,d+of],y(end)*[1 1],'b');
    
    p(2)=plot(x2,y3,'r--','Displayname','time average');
    p(7)=plot(xof,y3(end)*yof,'r--');
    
    
    p(3)=plot(x2,y2,'r.','Displayname','full/empty');
    p(8)=plot(xof,y2(end)*yof,'r.');
    p(4)=plot(x2,ones(size(x2)),'r.');
    p(9)=plot(xof,yof,'r.');
    
    xticks([0 X 1])
    xticklabels({'0','x','1'})
    yticks([0 0.5 1])
    xlh=xlabel('beam path');
    pub.BFxlab(xlh,T.F)
    ylh=ylabel('beam intensity');
    pub.BFylab(ylh,T.F)
    q(1)=quiver( 1.45,y3(end)+0.1,0,-0.1,0,'k','ShowArrowHead','on' )   ;
    q(2)=quiver( 1.45,y(end)-0.1,0,0.1,0,'k','ShowArrowHead','on' )   ;
    q(1).MaxHeadSize=10;
    q(2).MaxHeadSize=10;
    text(1.3,y3(end)+0.15,'bias')
    set(gca,'Layer','top')
   
    
    %text(0, 1.1,{'attenuating ','medium'})
    annotation('textarrow',[0.4 0.45],[0.35 0.4],'String',{'attenuating ','medium'})

    grid on
    xlim([-of,d+of])
    ylim([0,1.2])
    lh=legend(p(1:3),'location','southeast');
    pub.BFlegend(lh,T.F)
    
    fname='dnamicbias1';
    fpath=sprintf('%s%s%s',T.F.saveto,'technical\',fname);
    savefig(fh,fpath)
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpdf')
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpng')
%% 
    fh=figure(9876) ;
    clf;
    
      pub.BFfigure2(fh,T.F,1)
      fh.Position(1)=10;
    ax=gca;
    pub.BFaxis(ax,T.F)
    hold on
    
    of=0.5;
    n=100;
    d=2;
    x=linspace(0,1,100);
    mu=0.058;
    y=((x.*exp(-mu*d)+(1-x)))./exp(-mu*x*d)-1;
    
   
    p(1)=plot(x,y,'b','Displayname','bias');
    %p(2)=plot([-of,0],0*[1,1],'b');
    %p(3)=plot([d,d+of],y(end)*[1 1],'b');
    
    
    %xticks([0 X 1])
    %xticklabels({'0','x','1'})
    %yticks([0 0.5 1])
    ylh=ylabel('dynamic bias');
    xlh=xlabel('fill fraction x');
    pub.BFxlab(xlh,T.F)
    pub.BFylab(ylh,T.F)
    
    grid on
    xlim([0,1])
    %ylim([0,1.2])
    lh=legend(p(1));
    pub.BFlegend(lh,T.F)
    
    fname='dnamicbias2';
    fpath=sprintf('%s%s%s',T.F.saveto,'technical\',fname);
    savefig(fh,fpath)
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpdf')
    print(sprintf('%s%s%s',T.F.saveto,'technical\',fname),'-dpng')
    


    
end

function TomExample(T,F,mtomo)
    % the tomography vizualisation for presentation audicen
    % 1) make a phantom. use 269x269 in case we want to use real data later
    ph=zeros(269);
    phx=[120,180];
    phy=[150,50]
    %ph(200:210,50:80)=1;
    [x,y]=ndgrid(1:269,1:269);
    ph((x-phx(1)).^2+(y-phy(1)).^2<6^2)=0.5;
    ph((x-phx(2)).^2+(y-phy(2)).^2<20^2)=1;
    
    
    
    ph=ph/sum(ph(:)); %normalizing
    fid=4784;
    fh=figure(fid);    clf;
    
    xoff=0.05;
    p1x=1;
    p12x=0.1;
    p2x=0.5;
    p23x=0.1;
    p3x=360/269;
    
    yoff=0.05;
    p1y=1;
    yoff2=0.2
    
    
    fig_x_pos=[0 xoff p1x p12x p2x p23x p3x xoff];
    figxsum=sum(fig_x_pos);
    fig_x_pos=fig_x_pos/figxsum;
    fig_y_pos=[yoff p1y yoff2];
    figysum=sum(fig_y_pos);
    fig_y_pos=fig_y_pos/figysum;
    fh.Units='pixels'
    fscale=269; % pixels  and 80% coverage
    fh.Position=[0 0 figxsum figysum]*fscale;
    
    pos1=[(xoff)/figxsum yoff/figysum p1x/figxsum p1y/figysum];
    pos2=[(xoff+p1x+p12x)/figxsum yoff/figysum p2x/figxsum p1y/figysum];
    pos3=[(xoff+p1x+p12x+p2x+p23x)/figxsum yoff/figysum p3x/figxsum p1y/figysum];
    
    subplot('Position',pos1)
    imshow(ph,[])
    ang=0:359;
    pro=radon(ph,ang);
    pro2=pro(57+1:end-57,:); % cut edges
    proraw=zeros(size(pro2));
    maxpro=max(pro2(:));
    
    subplot('Position',pos2)
    xticks([])
    yticks([])
    
    subplot('Position',pos3)
    imshow(proraw,[]) 
   
    
    v = VideoWriter('V:\1phd\präsentation\PhD\forward');
    frate=20;
    v.FrameRate=frate;
        open(v)
    % rotate the phantom
    for a=1:length(ang)
        imrot=imrotate(ph,-ang(a),'bilinear','crop');
        subplot('Position',pos1)
        imshow(imrot',[]);set(gca,'YDir','normal');
        
        subplot('Position',pos2)
        plot(maxpro-pro2(:,a),1:269);
        xticks([])
        yticks([])
        xlim([0 maxpro])
        ylim([0 269])
        th=title(sprintf('%d°',a));
        th.FontSize=16;
        
        proraw(:,a)=pro2(:,a);
        subplot('Position',pos3)
        imshow(proraw,[]);set(gca,'YDir','normal');
 
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame)
    end

    close(v)
    
    %% the back projection.
    fid=5134;
    fh=figure(fid);    clf;
    
    xoff=0.05;
    p1x=p3x;
    p12x=0.1;
    p2x=1;
    p23x=0.1;
    p3x=1;
    
    yoff=0.05;
    p1y=1;
    yoff2=0.2
    
    
    fig_x_pos=[0 xoff p1x p12x p2x p23x p3x xoff];
    figxsum=sum(fig_x_pos);
    fig_x_pos=fig_x_pos/figxsum;
    fig_y_pos=[yoff p1y yoff2];
    figysum=sum(fig_y_pos);
    fig_y_pos=fig_y_pos/figysum;
    fh.Units='pixels'
    fscale=269; % pixels  and 80% coverage
    fh.Position=[0 0 figxsum figysum]*fscale;
    
    pos1=[(xoff)/figxsum yoff/figysum p1x/figxsum p1y/figysum];
    pos2=[(xoff+p1x+p12x)/figxsum yoff/figysum p2x/figxsum p1y/figysum];
    pos3=[(xoff+p1x+p12x+p2x+p23x)/figxsum yoff/figysum p3x/figxsum p1y/figysum];
    
    rec=zeros(269,269,length(ang));
    proj=zeros(269,269,length(ang));
    
    for a=1:length(ang) %compute recons
        a
        rec(:,:,a)=iradon(pro(:,1:a),ang(1:a),'spline','hann',269);
        proj(:,:,a)=iradon(pro(:,a),ang(a),'spline','hann',269);
    end
    rec(rec<0)=0;
    proj(proj<0)=0;
    
    v = VideoWriter('V:\1phd\präsentation\PhD\backward');
    v.FrameRate=frate;
    open(v)
    for a=1:length(ang) %compute recons
        subplot('Position',pos1);
        cla
        imshow(pro2,[]);set(gca,'YDir','normal');
        hold on
        plot([a,a],[0,269],'r')
        
        subplot('Position',pos2);
        imshow(imrotate(proj(:,:,a),ang(a),'bilinear','crop')',[]);set(gca,'YDir','normal');
        th=title(sprintf('%d°',a));
        th.FontSize=16;
        
        subplot('Position',pos3);
        imshow(rec(:,:,a)',[]);set(gca,'YDir','normal');
        
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame)
    end
    close(v)
    
end


function BFfigure(fh,F)
    % sets figure properties
    % fh = figure handle
    % F = Figure property struct (T.F)
    % use pub.BFfigure(h,T.F)
    set(fh,...
        'Units',F.Fig.Units,...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 2*F.Fig.FigW 2*F.Fig.FigH],...
        'PaperPosition',[0 0 2*F.Fig.FigW 2*F.Fig.FigH],...
        'PaperSize',2*F.Fig.PapS);
end

function BFfigure2(fh,F,x)
    % sets figure properties
    % fh = figure handle
    % F = Figure property struct (T.F)
    % x = zoom factor 12/16
    % use pub.BFfigure(h,T.F)
    set(fh,...
        'Units',F.Fig.Units,...
        'PaperPositionMode','manual',... %F.Fig.PPM,...
        'Position',[0 0 x*F.Fig.FigW x*F.Fig.FigH],...
        'PaperPosition',[0 0 x*F.Fig.FigW x*F.Fig.FigH],...
        'PaperSize',x*F.Fig.PapS);
end

function BFLine(p,F)
    % sets Line properties
    % p = line handle
    % F = Figure property struct (T.F)
    % use pub.BFLine(p,T.F)
    for i=1:length(p)
    set(p(i),...
        'Linewidth',F.L.LW);
    end
end

function BFxlab(xlh,F)
    % sets xlabel properties
    % xlh = xlabel handle
    % F = Figure property struct (T.F)
    % use pub.BFxlab(xlh,T.F)
    set(xlh,...
        'FontUnits',F.xl.FU,...
        'FontWeight',F.xl.FW,...
        'FontSize',F.xl.FS,...
        'FontName',F.xl.FN)
         %'Iterpreter',F.xl.Interp);
end

function BFylab(ylh,F)
    % sets ylabel properties
    % ylh = ylabel handle
    % F = Figure property struct (T.F)
    % use pub.BFylab(ylh,T.F)
    set(ylh,...
        'FontUnits',F.yl.FU,...
        'FontWeight',F.yl.FW,...
        'FontSize',F.yl.FS,...
        'FontName',F.yl.FN)
         %'Iterpreter',F.yl.Interp);
end

function BFaxis(ax,F)
    % sets axis properties
    % ax = axis handle
    % F = Figure property struct (T.F)
    % use pub.BFaxis(axlh,T.F)
    set(ax,...
        'FontUnits',F.ax.FU,...
        'FontWeight',F.ax.FW,...
        'FontSize',F.ax.FS,...
        'FontName',F.ax.FN,...
        'Units',F.ax.U);
        %'Iterpreter',F.ax.Interp
        
end

function BFtitle(th,F)
    % sets title properties
    % th = title handle
    % F = Figure property struct (T.F)
    % use pub.BFtitle(th,T.F)
    set(th,...
        'FontUnits',F.ti.FU,...
        'FontSize',F.ti.FS,...
        'FontName',F.ti.FN,...
        'FontWeight',F.ti.FW)
end

function BFlegend(lh,F)
    % sets legend properties
    % xlh = legend handle
    % F = Figure property struct (T.F)
    % use pub.BFlegend(lh,T.F)
    set(lh,...
        'FontUnits',F.L.FU,...
        'FontWeight',F.L.FW,...
        'FontSize',F.L.FS,...
        'FontName',F.L.FN)
        %'Iterpreter',F.L.Interp);
end


% pub.BFtitle(th,T.F)
% pub.BFfigure(fh,T.F)
% pub.BFfigure2(fh,T.F,x)
% pub.BFlegend(lh,T.F)
% pub.BFaxis(ax,T.F)
% pub.BFylab(ylh,T.F)
% pub.BFxlab(xlh,T.F)
% pub.BFLine(p,T.F)

end %static
end %class
















