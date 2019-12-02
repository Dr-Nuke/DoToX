function TomVid()
    % the tomography vizualisation for presentation audicen
    % 1) make a phantom. use 269x269 in case we want to use real data later
    tomsize=269;
    ph=zeros(tomsize);
    
    % add circles
    phx=[120,180];  % x-centers
    phy=[150,50];    % y-centers
    radii=[6 20];   % radii
    [x,y]=ndgrid(1:tomsize,1:tomsize);
    ph((x-phx(1)).^2+(y-phy(1)).^2<radii(1)^2)=0.5;
    ph((x-phx(2)).^2+(y-phy(2)).^2<radii(2)^2)=1;
    
    %add a slab
    %ph(200:210,50:80)=1;
    
    ph=ph/sum(ph(:)); %normalizing
    fid=4784;clf;
    fh=figure(fid);    clf;
    
    % "sophisticated" subplot-arrangment in units of full axis
    xoff=0.05;
    p1x=1;
    p12x=0.1;
    p2x=0.5;
    p23x=0.1;
    p3x=360/tomsize;
    
    yoff=0.05;
    p1y=1;
    yoff2=0.2;    
    
    fig_x_pos=[0 xoff p1x p12x p2x p23x p3x xoff];
    figxsum=sum(fig_x_pos);
    fig_x_pos=fig_x_pos/figxsum;
    fig_y_pos=[yoff p1y yoff2];
    figysum=sum(fig_y_pos);
    fig_y_pos=fig_y_pos/figysum;
    fh.Units='pixels'
    fscale=tomsize;
    fh.Position=[0 0 figxsum figysum]*fscale;
    
    % subplot positions
    pos1=[(xoff)/figxsum yoff/figysum p1x/figxsum p1y/figysum];
    pos2=[(xoff+p1x+p12x)/figxsum yoff/figysum p2x/figxsum p1y/figysum];
    pos3=[(xoff+p1x+p12x+p2x+p23x)/figxsum yoff/figysum p3x/figxsum p1y/figysum];
    
    %subplot('Position',pos1)
    %imshow(ph,[])
    ang=0:359;
    pro=radon(ph,ang);
    pro2=pro(57+1:end-57,:); % cut edges % hardcode warning
    proraw=zeros(size(pro2));
    maxpro=max(pro2(:));
    
    subplot('Position',pos2)
    xticks([])
    yticks([])
    
    subplot('Position',pos3)
    %imshow(proraw,[]) 
   
    
    v = VideoWriter('V:\1phd\präsentation\PhD\forward');
    frate=20; % increase for faster video
    v.FrameRate=frate;
        open(v)
    % xray line x coordinates
    linex=[0,tomsize];   
    nlines=6; % number of x-rays
    % forward projection    
    for a=1:length(ang)
        imrot=imrotate(ph,-ang(a),'bilinear','crop');
        subplot('Position',pos1)
        cla
        imshow(imrot',[]);set(gca,'YDir','normal');
        hold on
        
        % add xray lines
        
        for i=1:nlines
            liney=(i)/(nlines+1)*[1 1]*tomsize;
            plot(linex,liney,'y--')
            
            
            
        end
        
        subplot('Position',pos2)
        plot(maxpro-pro2(:,a),1:tomsize,'linewidth',2);
        xticks([])
        yticks([])
        xlim([0 maxpro])
        ylim([0 tomsize])
        th=title(sprintf('%d°',a));
        th.FontSize=16;
        
        proraw(:,a)=pro2(:,a);
        subplot('Position',pos3)
        imshow(proraw,[]);set(gca,'YDir','normal');
        hold on
        plot([a,a],[0,tomsize],'r','linewidth',2)
 
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
    
    
    fig_x_pos=[0 xoff p1x p12x p2x xoff];
    figxsum=sum(fig_x_pos);
    fig_x_pos=fig_x_pos/figxsum;
    fig_y_pos=[yoff p1y yoff2];
    figysum=sum(fig_y_pos);
    fig_y_pos=fig_y_pos/figysum;
    fh.Units='pixels'
    fscale=tomsize; % pixels  and 80% coverage
    fh.Position=[0 0 figxsum figysum]*fscale;
    
    pos1=[(xoff)/figxsum yoff/figysum p1x/figxsum p1y/figysum];
    pos2=[(xoff+p1x+p12x)/figxsum yoff/figysum p2x/figxsum p1y/figysum];
    %pos3=[(xoff+p1x+p12x+p2x+p23x)/figxsum yoff/figysum p3x/figxsum p1y/figysum];
    
    rec=zeros(tomsize,tomsize,length(ang));
    proj=zeros(tomsize,tomsize,length(ang));
    
    for a=1:length(ang) %compute recons
        a
        rec(:,:,a)=iradon(pro(:,1:a),ang(1:a),'spline','hann',tomsize);
        proj(:,:,a)=iradon(pro(:,a),ang(a),'spline','hann',tomsize);
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
        plot([a,a],[0,tomsize],'r','linewidth',2)
        
%         subplot('Position',pos2);
%         imshow(imrotate(proj(:,:,a),ang(a),'bilinear','crop')',[]);set(gca,'YDir','normal');
%         th=title(sprintf('%d°',a));
        
        
        subplot('Position',pos2);
        imshow(rec(:,:,a)',[]);set(gca,'YDir','normal');
        th=title(sprintf('%d°',a));
        th.FontSize=16;
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame)
    end
    close(v)
    
end