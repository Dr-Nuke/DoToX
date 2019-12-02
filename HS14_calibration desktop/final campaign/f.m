classdef f
    
properties
    T
end
%% function collection file for the final campaign Dec. 201+
methods(Static) %evil function scope hack

function [filenames,nframes,reps]=GenerateFileList(T,fid)
    % generates the files list in the folder, along with some parameters
    reps=zeros(T.d.ncas,T.d.nrep); % indicator of available files
    oldPath=cd(T.d.DataPath);
    list=dir('*tom_*.seq');
    cd(oldPath);
    filenames=cell(T.d.ncas,T.d.nrep);
    nframes=zeros(T.d.ncas,T.d.nrep);
    for i=1:length(list) % iterate mesurements
        fname=list(i).name;
        cas=str2num(fname(14:15));
        rep=str2num(fname(17:18));
        filenames{cas,rep}=list(i).name;
        nframes(cas,rep)=(list(i).bytes-T.d.header)...
            /T.d.GreyBytes/prod(T.d.imsize);
        if mod(nframes(cas,rep),1) ~= 0
            error('file size is not coherent with pixel sizes and header');
        end
        T.d.reps(cas,rep,1)=1;
    end
    
    h=figure(fid);clf
    d=nframes;
    d(d~=0)=d(d~=0)-min(d(d~=0));
    imagesc(d');colorbar()
    set(gca,'YDir','normal')
    ax=gca;
    title( f.FigTName2('difference of tomo lenghts (in frames)',0,0))
    ylabel('repetition')
    xlabel('case/ betriebspunkt')
    zlabel('additional frames')

    set(ax,'FontSize',14)
    tightInset = get(ax, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - 2.6*tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(ax, 'Position', position);
    set(ax,'units','centimeters')
    pos = get(ax,'Position');
    ti = get(ax,'TightInset');
    
    set(h, 'PaperUnits','centimeters');
    set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    fname=f.FigFileName2('Lengths differeces of raw tomos',0,0);

    print(h,fname,'-dpdf')
    print(h,fname,'-dpng')
    savefig(h,fname)
    
end

function [d]=ReadXrayTomo(T)
    % returns the xray raw data of a tomo
    % fpath = absolute path of the file to be read in
    % header = header size, usually 2048 for viva .seq videos
    % imsize = [x y] dimension, usually [640 1024] or [320 512] for 85Hz
    % bytes = color depth in bytes, usually 2

    fpath=strcat(T.d.DataPath,T.d.List{T.cas,T.rep});
    %fprintf('%s reading %s... ',f.CommLineStart(T.cas,T.rep),T.d.List{T.cas,T.rep})
    % open the file
    fid=fopen(fpath);

    % dispose of header
    fread(fid,T.d.header);%the header size migth need to be adjusted depending on image settings

    %read in data
    d=(fread(fid,T.d.imsize(1)*T.d.imsize(2)*T.d.nframes(T.cas,T.rep),'uint16=>single')); %the image
    fclose(fid); %close file handle

    % rotate such that the representation is correct
    d=permute(reshape(d,[T.d.imsize(2) T.d.imsize(1) T.d.nframes(T.cas,T.rep)]),[2,1,3]);
    %fprintf('done.')
end

function Dosemask=MakeDoseMask(T)
    % creates the dose mask
    Dosemask=zeros(T.d.imsize); %make the dose mask
    for i=1:size(T.dose.windows,2)
        Dosemask(T.dose.windows(1,i):T.dose.windows(2,i),...
               T.dose.windows(3,i):T.dose.windows(4,i))=1;
    end
end

function mask=MakeMask(windows,imsize)
    % creates the dose mask
    mask=zeros(imsize); %make the dose mask
    for i=1:size(windows,2)
        mask(windows(1,i):windows(2,i),...
               windows(3,i):windows(4,i))=1;
    end
    mask=logical(mask);
end

function d=SinoPadder(d,maxframes)
    % pads the frame dimension with images until the
    % total frame number is maxfraes
    if maxframes < size(d,3)
        error('maxframes is smaller than number of existing frames')
    end
    a=d(:,:,end);
    d=cat(3,d,repmat(a,[1,1,maxframes-size(d,3)]));
end

function []= f_BoFig2PDF(h,fname)
    % prints a figure to .fig,.pdf and .png
    % after resizing it
    % better use pub.PubPlot(gcf,gca,fname) nowadays
    set(h,'PaperOrientation','landscape');
    set(h,'PaperPosition', [1 1 28 19]);
    savefig(h,fname)
    print(h, '-dpdf',fname);
    set(h,'PaperOrientation','portrait');
    print(h, '-dpng',fname);
end

function []= f_BoFig3PDF(h,path,fname)
    % prints a figure to .fig,.pdf and .png
    % after resizing it
    % better use pub.PubPlot(gcf,gca,fname) nowadays
    set(h,'PaperOrientation','landscape');
    set(h,'PaperPosition', [1 1 28 19]);
    savefig(h,strcat(path,fname))
    %print(h, '-dpdf',strcat(path,fname));
    set(h,'PaperOrientation','portrait');
    print(h, '-dpng',strcat(path,fname));
end

function [ output_args ] = f_BoCount(c,f,r,d)
    % puts counter display values in order
    % f_BoCount(c,f,r,d)
    % c = counter
    % f = only every f-th number shown
    % r = output line width in numbers
    % d = digits per number

    if mod(c,f)==0
        string=strcat('%',num2str(d),'.d');
        if mod(c,r*f)==0
            fprintf(strcat(string,'\n'),c)

        else
            fprintf(string,c)
        end
    end

end

function figure01(im,cas,rep)
    h=figure(1);clf;
    subplot(2,1,2)
    imshow(squeeze(im),[]);
    fname=f.FigFileName('raw sinogram',1,T.cas,T.rep)
    f.f_BoFig2PDF(h,fname)
end

function BoSave(d,basename,T)
    fname=sprintf('%s%s%02d.mat',T.d.DataPath,basename,T.i);
    fprintf('saving %s...',fname)
    savefast(fname,'d');
    fprintf('done. \n')
end

function BoSave2(array,fname,fpath)
    fprintf('saving %s...',fname)
    savefast(strcat(fpath,fname,'.mat'),'array');
    fprintf('done. \n')
end

function fpath=BoSave3(d,basename,T)
    % saves according to cas & rep
    fname=sprintf('%s%02d_%02d.mat',basename,T.cas,T.rep);
    fpath=sprintf('%s%s',T.d.DataPath,fname);
    fprintf('%s saving %s...',f.CommLineStart(T.cas,T.rep),fname)
    savefast(fpath,'d');
    fprintf('done. \n')
end

function fpath=BoSave4(d,basename,T)
    %saves according to only cas
    fname=sprintf('%s%02d.mat',basename,T.cas);
    fpath=sprintf('%s%s',T.d.DataPath,fname);
    fprintf('%s saving %s...',f.CommLineStart2(T.cas),fname)
    savefast(fpath,'d');
    fprintf('done. \n')
end

function BoSaveSino(d,T)
    %save function specifically for the data fom 2017 12 13
    jmax=[10,10,9,10];
    for i=1:4
        for j=1:jmax
            fname=fprintf('case_%02d_%02d.mat',i,j);
            BoSave2(d,T.d.DataPath,fname)
        end
    end
end

function idiff=imdiff(im1,im2,x,y)
    a=abs(squeeze(im1-im2));
    idiff=sum(a(:));
end

function [shift,dif,shiftguess]=SinoMatch(ref,sino,T)
    %disp(sprintf('%02d %02d matching',T.cas, T.rep))
    % 1) do the rough finding via cross correlation
    %guessF=f.SinoMatchXCOV(ref,sino,T.Match.GuessLine);
    shiftguess=dftregistration(fft2(ref),fft2(sino),1000);

    % 2) do the fine matching iteratively
    %[shift,dif]=f.SinoMatchIterate(ref,sino,0,guessF,T);
    [shift,dif]=f.SinoMatchIterate(ref,sino,shiftguess(1),shiftguess(2),T);
    
    

end

function T=FindShiftDiffArea(T)
    %% manually find the shift diff area
    % just plots the mask
    h=figure(106);clf

    
    DiffAreaMask=zeros(size(T.sino.Raw(:,:,1)));
    imshow(squeeze(T.sino.Raw(:,:,1)),[])
    DiffAreaMask([10:190,460:630],50:1450)=1;
    imshow(0.5*DiffAreaMask+T.sino.Raw(:,:,1),[])
    T.Match.MaskSino=DiffAreaMask;
    fname=sprintf('shift area plot');
    f.f_BoFig2PDF(h,fname)
end

function shift=SinoMatchXCOV(s1,s2,pixrow)
    % the one row shift check
    [c,lags]=xcov(s1(pixrow,:),s2(pixrow,:),50,'coeff');
    if any(isnan(c)) %debug
        disp(sprintf('%4d contains NaN',i))
        shift=NaN;
    else % normal case
        %disp(sprintf('%4d does not contain NaN',i))
        [~,ind]=max(c);
        shift=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
            (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
    end
end

function shift=MaxXCOV(s1,s2)
    % general maximum xcov finder
    [c,lags]=xcov(s1,s2,'coeff');
    if any(isnan(c)) %debug
        disp(sprintf('%4d contains NaN',i))
        shift=NaN;
    else % normal case
        %disp(sprintf('%4d does not contain NaN',i))
        [~,ind]=max(c);
        shift=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
            (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
    end
end

function test_SinoMatchXCOV(s1,s2,pixrow)
    % call with
    % load('T')
    % f.test_SinoMatchXCOV(T.sino.Raw(:,:,1),T.sino.Raw(:,:,11),200)

    h=figure(101);clf;
    ax1=subplot(1,3,1);
    imshowpair(s1',s2');set(gca,'YDir','normal')
    line([pixrow,pixrow],[1,size(s1,2)])

    idiff=f.imdiff(s1,s2);
    shift=f.SinoMatchXCOV(s1,s2,pixrow);
    title(sprintf('diff= %.f',idiff));


    figure(102);clf;
    plot(s1(pixrow,:),'DisplayName','Img 1')
    hold on
    plot(s2(pixrow,:),'DisplayName','Img 2')
    title(sprintf('shift %f',shift))

    s21=f.fracshift(s2,[0,shift]);
    plot(s21(pixrow,:),'DisplayName','shifted Img 2')
    legend()

    figure(101);
    ax2=subplot(1,3,2);
    imshowpair(s1',s21');set(gca,'YDir','normal')
    line([pixrow,pixrow],[1,size(s1,2)])

    idiff=f.imdiff(s1,s21);
    title(sprintf('diff= %.f',idiff));

    ax3=subplot(1,3,3);
    imshow(s1'-s21',[]);set(gca,'YDir','normal')
    line([pixrow,pixrow],[1,size(s1,2)])
    title(sprintf('absolute difference'));

    linkaxes([ax1,ax2,ax3],'xy')
    fname='sino shift comparison';
    f.f_BoFig2PDF(h,fname)

end

function test_SinoMatchXCOV2(T)
    % call with
    % f.test_SinoMatchXCOV2(T)
    % x-dependecy
    pixrows=[1:640];
    h1=figure(105);clf;
    ax=gca;
    hold on
    grid on
    s1=T.sino.Raw(:,:,1);
    col=hsv(39);
    for j=2:39 %iterate data
        s2=T.sino.Raw(:,:,j);
        for i=pixrows %iterate rows
            shift(j,i)=f.SinoMatchXCOV(s1,s2,pixrows(i));
        end
        plot(shift(j,:),'color',col(j,:),'displayname',sprintf('dataset %d',j))
    end
    xlim([1,640])
    xlabel('x-pixel')
    ylabel('pixel shift')
    title(sprintf('dependency of pixel shift with respect to dataset 1\n on x-pixel rows '))
    
    fname=sprintf('pixel shift by XCOV');
    f.f_BoFig2PDF(h1,fname)
end

function diff2d=ShiftDiff(ref,sino,DiffMask)
    im=abs(ref-sino).*DiffMask
    
end

function [shift,dif]=SinoMatchIterate(ref,sino,guessX,guessF,T)
    % sino is being matched to ref
    
    %% iterate 2d shift grid
    yr=T.Match.YRange;
    sxb=linspace(-yr,yr,T.Match.YSteps);  % basis iteration grd
    sfb=linspace(-yr,yr,T.Match.YSteps);  % i is i-dimenion, f is frames-dimension
    sx=sxb+guessX;                 % actual iteration grid
    sf=sfb+guessF;
    fid=100;
    h=figure(fid);
    iterations=3;
    for j=1:iterations % number of iterations
        diff2d=zeros(length(sx),length(sf)); % preallocate map
        for k=1:length(sx) % x dimension
            for l=1:length(sf) % frame dimension
                m=l+(k-1)*length(sx);
                %f.f_BoCount(m,20,10,5)
                im=f.fraccircshift(sino,[sx(k),sf(l)]);
                
                % quadratic sum
                dif=(ref-im).^2.*T.Match.MaskSino;
                diff2d(k,l)=sum(dif(:));
            end
        end
        %fprintf('\n')
        % optimize for region
        [Min,I] = min(diff2d(:)); % find shift & index of minimum difference
        [I_row, I_col] = ind2sub(size(diff2d),I); % convert index to 2d index

        % plot
        ax(j)=subplot(1,iterations,j);
        imagesc(diff2d','XData', sx, 'YData', sf);set(gca,'YDir','normal');
        colorbar; colormap gray
        xlabel('x')
        if j==1
            ylabel('y')
        end

        % make title for plot
        titlestring=sprintf('iter. %d, min dev. %.1f \n shift: %4.2f %4.2f',...
            j,diff2d(I_row,I_col),sx(I_row),sf(I_col));
        if j==2
            ts=f.FigTName('',fid,T.cas,T.rep);
            titlestring=sprintf('%s \n %s',ts,titlestring);
        end
        title(titlestring)
        axis equal
        axis tight
        shift=[sx(I_row),sf(I_col)];
        dif=diff2d(I_row,I_col);
        hold on
        plot(sx(I_row),sf(I_col),'xr')
        
        % zoom in the shifts
        sxb=sxb/10;
        sfb=sfb/10;
        sx=sxb+sx(I_row);
        sf=sfb+sf(I_col);
    
    end
    fname=f.FigFileName('sino shift match',fid,T.cas,T.rep);
    f.f_BoFig2PDF(h,fname)
    
    cminmax=0.1*[-1,1];
    fid=110;
    h=figure(fid);clf;
    subplot(2,1,1)
    im=f.fracshift(sino,shift);
    imshow((ref-im).*(0.5+0.5*T.Match.MaskSino),[]);
    dif=abs(ref-im).*T.Match.MaskSino;
    dif=sum(dif(:));
    
    ts=sprintf('difference after shifted (%.2f,%.2f), diff %.f',shift(1),shift(2),dif);
    title(f.FigTName(ts,fid,T.cas,T.rep))
    colorbar
    caxis(cminmax)
    subplot(2,1,2)
    im2=f.fracshift(sino,[0,shift(2)]);
    imshow((ref-im2).*(0.5+0.5*T.Match.MaskSino),[]);
    dif2=abs(ref-im2).*T.Match.MaskSino;
    dif2=sum(dif2(:));
    title(sprintf('without x-shift, dif %.f',dif2))
    colorbar
    caxis(cminmax)
    fname=f.FigFileName('sinogram match',fid,T.cas,T.rep);
    f.f_BoFig2PDF(h,fname)
    
    %M.shift(i,:)=[]
end

function [shift,mindif]=SinoMatchY(T)
    % special shift to correct for an y-offset in HS17 Xray Tomo
    % run test-yshift to check upon this routine

    fid=120;
    h=figure(fid);clf;
    sb=linspace(-T.Match.YRange,T.Match.YRange,T.Match.YSteps);
    s=sb+T.Match.YGuess;
    iterations=3;
    for j=1:iterations % number of iterations
        dif2=zeros(length(s),1); % preallocate map
        for l=1:length(s) % frame dimension
            im=f.fraccircshift(double(squeeze(T.Match.YPlanes(:,:,T.cas,T.rep))),[s(l),0]);
            % dif=abs(double(squeeze(T.Match.YPlanes(:,:,1,1)))-im).*(T.Match.MaskY);
            dif=((double(squeeze(T.Match.YPlanes(:,:,1,1)))-im).*(T.Match.MaskY)).^2;
            dif2(l)=sum(dif(:));
        end
        [Min,I] = min(dif2); % find shift & index of minimum difference
        I_row = ind2sub(size(dif2),I); % convert index to 2d index
        shift=s(I_row);
        mindif=dif2(I_row);
        
        %plot
        subplot(1,3,j)
        plot(s,dif2,'-x')
        grid on
        xlabel('y-shift')
        if j==1
            ylabel('difference')
        end
        % make title for plot
        titlestring=sprintf('iter. %d, min dev. %.1f \n shift: %4.2f ',...
            j,mindif,shift);
        if j==2
            titlestring=sprintf('case %d %d \n %s',T.cas,T.rep,titlestring);
        end
        title(titlestring)
        % zoom in the shifts
        sb=sb/10;
        s=sb+shift;
    end
    
    fname=f.FigFileName('y-shift correction',fid,T.cas,T.rep);
    f.f_BoFig2PDF(h,fname)
    
    
end

function sino=fracshift(sino,shift)
    int = floor(shift);     %integer portions of shiftsize
    fra = shift - int;      %fractional portions of shiftsize
    sino=(1-fra(1))*f.sinoshift(sino,[int(1),0])... % x-shift
        +fra(1)*f.sinoshift(sino,[int(1)+1,0]);
    sino=(1-fra(2))*f.sinoshift(sino,[0,int(2)])... % F-shift
        +fra(2)*f.sinoshift(sino,[0,int(2)+1]);
    
    
end % fractional 2d shift

function test_fracshift(im,shift)
    % tests the fracshift function
    h=figure(103);clf
    imax=9;
    ishift=linspace(-shift,shift,imax);
    for i=1:imax
        subplot(3,3,i)
        %imshowpair(im',f.fracshift(im,[-ishift(i),ishift(i)])');
        imagesc(f.fracshift(im,[-ishift(i),ishift(i)])')
        set(gca,'YDir','normal')
        title(sprintf('shift %.1f',ishift(i)))
        grid on
    end
    fname='test_fracshifter'
    f.f_BoFig2PDF(h,fname)
end

function sino=sinoshift(sino,shift)
    % shifts a sinogram by shift pixels and pads the created area with
    % nearest neigbor
    if all(mod(shift,1))==0; %check if it is integer shift
        % x-shift
        if shift(1)>0
            sino=cat(1,repmat(sino(1,:),[abs(shift(1)),1]),...
                sino(1:end-shift(1),:));
        elseif shift(1)<0
            sino=cat(1,sino(-shift(1)+1:end,:),...
                repmat(sino(end,:),[abs(shift(1)),1]));
        end        
        % frames-shift
        if shift(2)>0
            sino=cat(2,repmat(sino(:,1),[1,abs(shift(2))]),...
                sino(:,1:end-shift(2)));
        elseif shift(2)<0
            sino=cat(2,sino(:,-shift(2)+1:end),...
                repmat(sino(:,end),[1,abs(shift(2))]));
        end
    else
        error('shift is not integer')
    end
end % integer 2d shift

function B = fraccircshift(A,shiftsize)
    %fraccircshift expands circshift to fractional shifts values, using linear
    %interpolation. In contrast to other approaches to non-integer shifts
    %of matrices on the base of fft2 or interp2, the number of dimensions of A
    %is not limited, and the syntax of circshift applies. For integer elements
    %of shiftsize, fracircshift and circshift give the same results.
    
    int = floor(shiftsize);     %integer portions of shiftsize
    fra = shiftsize - int;      %fractional portions of shiftsize
    dim = numel(shiftsize);
    B = A;
    for n = 1:numel(shiftsize)  %The dimensions are treated one after another.
        intn = int(n);
        fran = fra(n);
        shift1 = zeros(dim,1);
        shift1(n) = intn;
        shift2 = zeros(dim,1);
        shift2(n) = intn+1;
        %Linear intepolation:
        B = (1-fran)*circshift(B,shift1) + fran*circshift(B,shift2);
    end
end

function B = gpufraccircshift(A,shiftsize)
    % same as fraccircshift, but supporting GPU execution
    
    int = floor(shiftsize);     %integer portions of shiftsize
    fra = shiftsize - int;      %fractional portions of shiftsize
    dim = numel(shiftsize);
    B = A;
    for n = 1:numel(shiftsize)  %The dimensions are treated one after another.
        intn = int(n);
        fran = fra(n);
        shift1 = zeros(dim,1,'gpuArray');
        shift1(n) = intn;
        shift2 = zeros(dim,1,'gpuArray');
        shift2(n) = intn+1;
        %Linear intepolation:
        B = (1-fran)*circshift(B,shift1) + fran*circshift(B,shift2);
    end
end

function test_sinoshift(sino,shift)
    % shows how sinoshift works
    h=figure(102);clf;
    x=[-shift(1),0,shift(1)];
    y=[-shift(2),0,shift(2)];
    for i=1:3
        for j=1:3
            k=j+(i-1)*3;
            disp([i,j,k])
            subplot(3,3,k)
            imshow(f.sinoshift(sino,[x(i),y(j)])');set(gca,'YDir','normal')
            if k==5 %(middle one)
                title('original')
            else
                title(sprintf('shift (%d, %d)',x(i),y(j)))
            end
        end
    end
    fname='test_sinoshifter';
    f.f_BoFig2PDF(h,fname)

end

function nframes=frames360(sino,T)
    % returns the number of frames in a sinogram that results
    % in the 360° cut
    
    % we pick a an x-row of pixels, and find the most similar at a later
    % point
    cas=T.cas;
    rep=T.rep;
    refrow=T.q360.startframe;
    g=T.q360.NFGuess;
    r=T.q360.NFGrange;
    s=size(sino);
    l=squeeze(sino(:,refrow)); % the x-line
    for i=1:s(2)
        diff(i)=sum(abs(l-squeeze(sino(:,i))));
    end
    ld=length(diff);
    fr=1:(refrow+g-r);%front range
    
    if ld>refrow+g+r
        br=(refrow+g+r):ld;
    else
        br=[];
    end
    diff2=diff;
    diff([fr,br])=max(diff);
    [mindif,mini]=min(diff);
    nframes=mini-refrow;
    
    fid=20;
    h=figure(fid);clf
    plot(diff2,'DisplayName','consecutive difference')
    hold on
    plot(diff,'DisplayName','guess-crop difference')
    ts=f.FigTName(sprintf('finding the 360° frames\n nframes=%d',nframes),fid,cas,rep);
    title(ts)
    xlabel('frame number')
    ylabel('difference')
    grid on
    legend('Location','south')
    fname=f.FigFileName('crop 360 degree',fid,cas,rep);
    f.f_BoFig2PDF(h,fname)    
    
    fid=30;
    h=figure(fid);clf;
    ax1=subplot(3,1,1);
    
    start=T.Raw.CropRange(cas,rep,1,2)+refrow; % image start row
    stop=refrow+nframes+T.Raw.CropRange(cas,rep,1,2)
    
    imshow(circshift(sino(:,start:stop-2),[0,50]),[])
    title(f.FigTName(sprintf('\n -1 frame'),fid,cas,rep))
    
    

    ax2=subplot(3,1,2);
    imshow(circshift(sino(:,start:stop-1),[0,50]),[]) 
    title('correct cut')
    ax3=subplot(3,1,3);
    imshow(circshift(sino(:,start:stop-0),[0,50]),[])    
    title('+1 frame')
    linkaxes([ax1,ax2,ax3])
    xlim([20,80])
    ylim([380,430])
    
    
    fname=f.FigFileName('crop 360 degree',fid,cas,rep);
    f.f_BoFig2PDF(h,fname)
end

function [fits,fitshift,centershift]=FindCentering(T,d)
    % finds the per plane shift required for centering
    % and fits a linear function over the channel height
    % T = T-struct
    % d = raw data block, such that d(:,xx,:) is a sinogram
    % fits = the parameters of the fit
    % fitsshift = the fitted shifts
    % centershift = the calculated shifts
    
    cas=T.cas;
    rep=T.rep;
    k180=T.q360.nFrames(1,1)/2;
    k180frac=rem(k180,1);
    
    centershift=zeros(1,T.Raw.BS(2));
    fprintf('%s NaNs in ',f.CommLineStart2(cas))
    
    for plane = 1:T.Raw.BS(2) % iterate planes
        % determine centering
        row180=(1-k180frac)*d(:,plane,floor(k180))+...
            k180frac*d(:,plane,ceil(k180));
        [c,lags]=xcov(d(:,plane,1),flipud(row180),50,'coeff');
        if any(isnan(c)) %debug
            fprintf('%d, ',plane)
            continue
        else
            %disp(sprintf('%4d does not contain NaN',j))
            %ind_nan(k)=1;
        end
        [~,ind]=max(c);
        try
        centershift(plane)=(lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
            (log(c(ind-1))-2*log(c(ind))+log(c(ind+1))))/2;
        catch
            centershift(plane)=nan;
        end
    end
    fprintf('\n')
    
    % decide which points are used for the fit
    fitbase=zeros(1,T.Raw.BS(2),'logical');
    if 1%cas > 1
        fitbase(T.Cen.CSFR{2})=1;
    else
        fitbase(T.Cen.CSFR{1})=1;
    end
    fitbase(isnan(centershift))=0;
    planes=1:T.Raw.BS(2);
    %fits=polyfit(T.Cen.CSFR{2},centershift(T.Cen.CSFR{2}),1);
    fits=polyfit(planes(fitbase),centershift(fitbase),1);
    fitshift=polyval(fits,1:T.Raw.BS(2));
    %T.c.centershift{i}=centershift;
    if T.proofs==1
        f.CheckCenterFit(d,fits,fitshift,centershift,fitbase,T);
    end
    
end

function CheckCenterFit(d,fits,fitshift,centershift,fitbase,T)
    planes=1:length(centershift);
    fid=130;
    h=figure(fid);clf;
    p=plot(planes(fitbase),centershift(fitbase),'Displayname', 'fit base',...
        'Linewidth',6);
    p.Color(4) = 0.5;
    hold on
    plot(centershift,'Displayname', 'per-plane shift')
    
    plot(fitshift,'Displayname','fitted shift')
    xlim([1,1024])
    ylim([min(fitshift)-5,max(fitshift)+5])
    %ylim([-10,5])
    title(sprintf('%s\n%s',f.FigFileName(...
        'centering per y-plane and linear fit',fid,T.cas,T.rep),...
        sprintf('shift from %.1f to %.1f, diff %.3f',...
        fitshift(1),fitshift(end),fitshift(end)-fitshift(1))))
    xlabel('y plane (pixel planes)')
    ylabel('centering shift in pixel')
    grid on

    legend()
    fname=f.FigFileName('centering shift',fid,T.cas,T.rep);
    f.f_BoFig2PDF(h,fname)
end

function [ sino ] = f3_center( sino,shift )
    % this function performs the sub pixel shift with a given sinogram
    inte=fix(shift);
    frac=shift-fix(shift);
    startsize=size(sino,1);

    if fix(shift)>0 %in case of positive shift value, shifting to the left
        sino=[sino;repmat(sino(end,:),abs(inte),1)]; %integer shift to center
        bs=size(sino); %blocksize
        
        if mod(size(sino,1),2) % if the number of pixels is odd
            [XI,YI]=ndgrid(1:bs(1),1:bs(2));
            sino=interpn(XI,YI,sino,XI+double(frac),YI,'cubic'); %subgrid shift
        else % if it is even, make it odd
            [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
            sino=interpn(XI,YI,[sino;sino(end,:)],XI+double(frac)-0.5,YI,'cubic'); %subgrid shift
        end %this is to get always odd sized projections, that preserves centeredness in iradon
        
    else %in case of negative shift value, shifting axis to the right
        sino=[repmat(sino(1,:),abs(inte),1);sino]; %integer shift to center
        bs=size(sino); %blocksize
        
        if mod(size(sino,1),2) % if the number of pixels is odd
            [XI,YI]=ndgrid(1:bs(1),1:bs(2));
            sino=interpn(XI,YI,sino,XI+double(frac),YI,'cubic'); %subgrid shift
        else % if it is even, make it odd
            [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
            sino=interpn(XI,YI,[sino;sino(end,:)],XI+double(frac)-0.5,YI,'cubic'); %subgrid shift
        end %this is to get always odd sized projections, that preserves centeredness in iradon
    end
    endsize=size(sino,1);
    d=(endsize-startsize)/2;
    if d~=0
        if mod(d,1)==0
            sino=sino(1+d:end-d,:);
        else
            error(sprintf('wrong padding'))
        end
        
    end
        
end

function [ sino4 ] = SinoCentForRec( sino,shift )
    % this function performs the sub pixel centering before reconstruction.
    % sino = sinogram to be centered. must have odd width
    % shift = shift distance
    
    if mod(size(sino,1),2)==0;
        error('sinogram has even width. needs to be odd')
    end
    sizebak=size(sino);
    inte=fix(shift);
    frac=shift-fix(shift);

    
    % pad sino on both sides such that it wont give nans when
    % interpolationg
    npad=ceil(abs(shift)); 
    xp=5; %extra pad
    sino2=[repmat(sino(1,:),npad+xp,1);sino;repmat(sino(end,:),npad+xp,1)];
    ss=size(sino2);
    
    [xx,yy]=ndgrid(1:ss(1),1:ss(2));
    sino3=interpn(xx,yy,sino2,xx+shift,yy,'cubic');
    
    sino4=sino3((npad+xp+1):(end-(npad+xp)),:);
   
end

function nframes=FindNframes(T,i)
    % returns the number of frames of .seq a video file
    fpath=strcat(T.d.DataPath,T.d.List(i).name);
    FileInfo=dir(fpath);
    FileSize=FileInfo.bytes;
    nframes=(FileSize-T.d.header)/(T.imsize(1)*T.imsize(2)*T.d.GreyBytes);
    if mod(nframes,1) ~= 0
        error('file size is not coherent with pixel sizes and header');
    end
end

function maxframes=FindMaxFrames(T)
    % this function does not work because the max frames number is
    % assessed only after sinograms of all data were made
    nframes=zeros(1,size(T.d.List,1));
    for i = 1:size(T.d.List,1)
        nframes(i)=f.FindNframes(T,i);
    end
    maxframes=max(nframes);
end

function y = loadSingleVariableMATFile(filename)
    % used to load a saved mat array into a spcifiable variable
    fprintf('loadind %s .....',filename)
    foo = load(filename);
    whichVariables = fieldnames(foo);
    if numel(whichVariables) == 1
        y = foo.(whichVariables{1});
        fprintf('done\n')
    else
        error('the file loaded contains another number of variables than 1. nothing was loaded.')
    end
end

function y = BoLoad(filename,T)
    [filepath,name,ext] = fileparts(filename);
    fprintf('%s loadind %s .....',f.CommLineStart2(T.cas),name)
    foo = load(filename);
    whichVariables = fieldnames(foo);
    if numel(whichVariables) == 1
        y = foo.(whichVariables{1});
        fprintf('done.\n')
    else
        error('the file loaded contains another number of variables than 1. nothing was loaded.')
    end
end

function flt=mediflt(im)
    flt=medfilt2(im, [3 3])
end

function [fitresult, gof] = FitPixel1(xring, yring, badimring)
    % firts the detector bad pixel things
%CREATEFIT(XRING,YRING,BADIMRING)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xring
%      Y Input : yring
%      Z Output: badimring
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Jan-2018 12:20:44


%% Fit: 'untitled fit 1'.
%[xData, yData, zData] = prepareSurfaceData( -yring, xring, badimring );
% bakup:
[xData, yData, zData] = prepareSurfaceData( xring, yring, badimring );
% Set up fittype and options.
ft = fittype( 'poly22' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'badimring vs. xring, yring', 'Location', 'NorthEast' );
% % Label axes
% xlabel xring
% ylabel yring
% zlabel badimring
% grid on
% view( -79.9, 11.6 );
end    

function im_filt=RemoveCameraSpots(im,cas,rep)
    %creates a mask that corrects the camera specific spots
    % manual spot findgin
    warning('off', 'curvefit:fit:equationBadlyConditioned')
    
    % 1) two dark spots x&y, box & center
    im_filt=im;
    dscxy = [294,405;... % center x
             473,567];  % center y
    mask=~mt.KernelGenBin(11,5,0);
    % visualize mask:
    % figure(126);clf;imagesc(mask)
    fid=80; 
    
    h=figure(fid);clf;
    for spot=1:size(dscxy,2)
        % create ringed data
        areax=(dscxy(1,spot)-5):(dscxy(1,spot)+5);
        areay=(dscxy(2,spot)-5):(dscxy(2,spot)+5);
        [x,y]=ndgrid(areax,areay);
        data=im(areax,areay);
        xr=x(mask);
        yr=y(mask);
        datar=double(data(mask));
        
        xn=x(~mask);
        yn=y(~mask);
        [fitresult{spot}, gof{spot}] = f.FitPixel1(xr, yr, datar);
        subplot(2,1,spot);
        p = plot( fitresult{spot}, [xr, yr], datar );
        set(p(1),'FaceAlpha',0.2)
        hold on
        datan=data(~mask);
        scatter3(xn,yn,datan,'MarkerFaceColor',[1,0,0]);
        im_filt(sub2ind(size(im), xn, yn))=feval(fitresult{spot},xn,yn);
        view( -79.9, 20 );
    end
    
  
    %legend( h, 'untitled fit 1', 'badimring vs. xring, yring', 'Location', 'NorthEast' );
    % Label axes
    xlabel x
    ylabel y
    zlabel counts
    
    grid on
    view( -79.9, 11.6 );
    warning('on', 'curvefit:fit:equationBadlyConditioned')
    subplot(2,1,1)
    legend('fit basis','faulty data','location','southwest')
    title(f.FigTName('camera spot fixes',fid,cas,rep))
    fname=f.FigFileName('camera spot fixes',fid,cas,rep);
    f.f_BoFig2PDF(h,fname)
    
end

function [Raw,bot,up,rad,yplane]=PreRead(T)
    toc1=toc;
    % reads in all data, and provides sinograms & radiograms
    % for manual parameter adjustment
    s=size(T.d.List);
    is=T.d.imsize;
    sy=T.Raw.sinoyplanes;
    Raw=cell(10,10);
    bot=cell(10,10);
    up=cell(10,10);
    rad=cell(10,10);
    yplane=cell(10,10);
    %preallocate
    maxcas=s(1);
    maxrep=s(2);
    
    for cas=1:maxcas
        for rep=1:maxrep
            if T.d.nframes(cas,rep) ~=0
                ind=rep+(cas-1)*maxrep;
                %disp([cas,rep,ind])
                d=f.ReadXrayTomo(T,cas,rep);
                
                % keep a copy of sinograms
                Raw{cas,rep}=squeeze(d(:,sy(1),:)); % corrected raw sinogram with angular marker
                bot{cas,rep}=squeeze(d(:,sy(2),:)); % bottom sinogram, above marker
                up{cas,rep}=squeeze(d(:,sy(3),:)); % high up sinogram
                rad{cas,rep}=squeeze(d(:,:,10)); %radiogram
                yplane{cas,rep}=squeeze(d(T.Raw.yshiftplane,:,:));
                %f.figure01(Raw{cas,rep},cas,rep); %make figure & save
            end
        end
    end
    toc2=f.NiceTime(toc1);
    fprintf('Elapsed time of PreRead(): %s \n',toc2)
end

function s=NiceTime(toc1)
    toc2=toc;
    % displays time difference in a nice way
    
    t=toc2-toc1;
    if t < 120
        s=sprintf('%.2f %s',t,'s');
    elseif t<7200
        s=sprintf('%.2f %s',t/60,'min');
    else
        s=sprintf('%.2f %s',t/3600,'h');
    end
    
end

function [croprange]=FindCrop(T,sino)
    % finds start and end of the raw rotation (>360°!)
    % logic: first and last time the minimum difference between two
    % consecutive frames is larger than a threshold
    croprange=zeros(2,3);
    a=sino(:,2:end)-sino(:,1:end-1);
    mi=min(a,[],1);
    
    % find the flanks 
    % start flank
    croprange(1,1)=find(abs(mi)>T.Raw.CropThrash,1,'first')-T.Raw.CropMargin;
    % stop plank; % make upper limit % that -4 i saw empirically: 4 frames
    % of zero at the end screwed up this algo
    approx=min(croprange(1,1)+T.Raw.CropMargin+T.Raw.RotCutLimit,T.d.nframes(T.cas,T.rep)-5);
    croprange(2,1)=find(abs(mi(1:approx))>T.Raw.CropThrash,1,'last')+T.Raw.CropMargin;
    
    % now, we might have found different long sinograms than the [1,1] case
    % lets fix that by assuming the difference is smaller than 
    % T.Raw.CropMargin, and remove a few pixels on both sides
    
    Len=croprange(2,1)-croprange(1,1)+1;
    Len_ref=T.Raw.CropRange(1,1,2,2)-...
             T.Raw.CropRange(1,1,1,2)+1;
    diff=Len_ref-Len;
    if mod(diff,1)~=0 % its not an integer!
        error('sino length diff is not integer! rc=(%.2f,%.2f), diff=%.2f',rc(1),rc(2),diff)
    end
    a=ceil(diff/2);
    b=floor(diff/2);
    croprange(:,2)=croprange(:,1)+[-a;b]; 

    % now, the croprange values could be before the first or after the last
    % frame.lets check. eventually, pad the data block accorfingly
    if and(croprange(1,2)<1,croprange(2,2)>T.d.nframes(T.cas,T.rep))
        error(sprintf('cant pad in both directions. cas %d rep %d'),T.cas,T.rep)
    end
    
    if croprange(1,2)<1 % need to pad ath the beginning
        % pad by the difference
        croprange(1,3)=1-croprange(1,2); % the padding length 
        %d=cat(3,repmat(d(:,:,1),[1,1,padLen]),d);
        %correct rc
        croprange(:,2)=croprange(:,2)+croprange(1,3);
    end
        
    if croprange(2,2)>T.d.nframes(T.cas,T.rep) % need to pad at the end
        % pad by the difference
        croprange(2,3)=croprange(2,2)-T.d.nframes(T.cas,T.rep);
        %d=cat(3,d,repmat(d(:,:,end),[1,1,padLen]));
        croprange(:,2)=croprange(:,2)-croprange(2,3);
    end
    
    if any(croprange(:,3)<0)
        error(sprintf('crop padding is negative. cas %d rep %d'),T.cas,T.rep)
    end
end

function CheckSino1(T,d,fid)
if T.proofs==1
  
        h=figure(fid);clf;
        % raw sinogram plot
        subplot(2,1,1)
        imshow(squeeze(d(:,T.Raw.sinoyplanes(1),:)),[]);
        ax=gca;
        hold on
        for i=1:2
            line(T.Raw.CropRange(T.cas,T.rep,i,2)*[1,1],...
                [1,T.d.imsize(2)],'color',[1,0,0])
        end
        tstr1='Raw sinogram of the angle gauge region';
        tstr2='with automated crop (red), cropped below';
        title(sprintf('%s \n %s',f.FigTName(tstr1,fid,T.cas,T.rep),tstr2));

    end 
end

function CheckSino2(T,d,fid)
     if T.proofs==1
        
        h=figure(fid);

        subplot(2,1,2);
        imshow(squeeze(d(:,T.Raw.sinoyplanes(1),:)),[]);
        ax=gca;
        
        subplot(2,1,1)

        line(T.q360.startframe*[1,1],[1,T.d.imsize(2)],'color',[0,1,0])
        line((T.q360.nFrames(T.cas,T.rep)+T.q360.startframe)*[1,1],...
            [1,T.d.imsize(2)],'color',[0,1,0])
        
        fname=sprintf('Fig01 %02d %02d raw sinogram',T.cas,T.rep);
        fname=f.FigFileName('raw sinogram',fid,T.cas,T.rep);

        f.f_BoFig2PDF(h,fname)
        
    end 
end

function [d,nframes]=SinoCrop(d,T)
    % acrop the datablock. pads if neccesseary
    
    if T.Raw.CropRange(T.cas,T.rep,1,3)>0
        d=cat(3,repmat(d(:,:,1),[1,1,T.Raw.CropRange(T.cas,T.rep,1,3)]),d);
    end
    
    if T.Raw.CropRange(T.cas,T.rep,2,3)>0
        d=cat(3,d,repmat(d(:,:,end),[1,1,T.Raw.CropRange(T.cas,T.rep,2,3)]));
    end
    
    d=d(:,:,T.Raw.CropRange(T.cas,T.rep,1,2):T.Raw.CropRange(T.cas,T.rep,2,2));
    nframes=size(d,3);
end

function [d_out,fi]=BeamAndDoseCorrection2(d,T)
    % corrects for beam & dose
    % d = data block
    % T = parameters struct
    %fprintf('%d %d correctiong for beam and dose...\n',T.cas,T.rep)
    mask=T.dose.mask; %re-label
    
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single'); % contains the correction factor

    % find the fits
    %fprintf('%d %d finding beam non-uniformity fits...\n',T.cas,T.rep)
    nframes=T.Raw.nframes(T.cas,T.rep); % parfor-hack
    parfor j=1:nframes

        %f.f_BoCount(j,200,10,5)
        im=squeeze(d(:,:,j));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fi{j}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
        unibeam(:,:,j)=feval(fi{j},x,y);         % make the beam

    end
   
    %fprintf('\n')

    % correct the non uniformity
    d_out=d./unibeam;

    % Dose Area checks
    if T.proofs==1 % check the areas
        fid=40;
        h=figure(fid);clf;
        set(gcf,'name','dose area mask check')
        subplot(1,3,1)
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        title('original')
        subplot(1,3,2)
        imshow(T.dose.mask',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d mask',T.cas,T.rep))
        subplot(1,3,3)
        imshow((single(squeeze(d(:,:,1))).*(0.5+T.dose.mask/2))',[]);set(gca,'YDir','normal')
        title(f.FigTName('combined',fid,T.cas,T.rep))
        title('combined')
        fname=f.FigFileName('dose mask',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)



        %check if nothing rotates in these areas
        fid=50;
        h=figure(fid);clf;
        dosemin=zeros(1,T.Raw.nframes(T.cas,T.rep));
        dosemean=zeros(1,T.Raw.nframes(T.cas,T.rep));
        for i=1:T.Raw.nframes(T.cas,T.rep)
            a=single(squeeze(d(:,:,i))).*T.dose.mask;
            dosemin(i)=min(a(a>0));
            dosemean(i)=mean(a(a>0));
        end

        subplot(2,1,1)
        plot(dosemean)
        title(f.FigTName('mean dose in dose areas per frame',fid,T.cas,T.rep))

        ylabel('mean dose')
        xlabel('frame number')
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        grid on

        subplot(2,1,2)
        plot(dosemin)
        ylabel('minimum dose image value')
        xlabel('frame number')
        title(f.FigTName(sprintf(...
            'if this has major dips, then some object \n rotates into the dose mask'),...
            fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)


    end

    %fit checks
    if T.proofs==1
        fid=60;
        h=figure(fid);clf;
        imax=1;
        ii=ceil(T.Raw.nframes(T.cas,T.rep)*rand(1,imax)); %check some random frames
        for i = 1:length(ii)

            im=squeeze(d(:,:,ii(i)));
            im(T.dose.mask==0)=[];
            rf =10; % reduction factor
            plot(fi{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
            hold on
            %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
            xlabel('x')
            ylabel('y')
            zlabel('counts')
            %pause(0.5)
        end
        view([20,45])
        title(f.FigTName(sprintf(...
            'nun-uniformity data \n according to a random frame'),...
            fid,T.cas,T.rep))
        fname=f.FigFileName('Beam fit',fid,T.cas,T.rep);    
        f.f_BoFig2PDF(h,fname)

        %check if the corrected values are ok
        fid=70;
        h=figure(fid);clf;
        subplot(2,3,4:5);
        hold on
        ndf=size(T.dose.windows,2); % number of dose fields

        doseprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        oldprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        dosepix=zeros(1,ndf);

        for i=1:ndf %iterate dose fields
            % mean of each field for each frame, for both old (im) and
            % corrected (imc) data
            doseprofile(:,i)=squeeze(mean(mean((d_out(...
                T.dose.windows(1,i):T.dose.windows(2,i),...
                T.dose.windows(3,i):T.dose.windows(4,i),:)))));
            oldprofile(:,i)=squeeze(mean(mean((d(...
                T.dose.windows(1,i):T.dose.windows(2,i),...
                T.dose.windows(3,i):T.dose.windows(4,i),:)))));

            % plot
            plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))

            % the number of dose mask pixels
            dosepix(i)=(T.dose.windows(2,i)-T.dose.windows(1,i)+1)...
                *(T.dose.windows(4,i)-T.dose.windows(3,i)+1);

            % the total dose, should be constantly equal 1
            doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
            oldprofile(:,ndf+1)=oldprofile(:,ndf+1)+dosepix(i)*oldprofile(:,i);
        end
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
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
        title(f.FigTName('mean per frame dose of raw footage',fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])

        subplot(2,3,[3,6]);
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        hold on
        for i=1:ndf
            r=rectangle('Position',[T.dose.windows(1,i),...
                T.dose.windows(3,i),...
                T.dose.windows(2,i)-T.dose.windows(1,i),...
                T.dose.windows(4,i)-T.dose.windows(3,i)],...
                'Facecolor','r');
        end
        title('dose areas')
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
        


    end
   
end

function [d_out,fdose,fdosemin]=DoseCorrect(d,T)
    % corrects only the dose, framewise
    % d_out = the dose corrected block
    % fdose = the per frame dose information
    % d = the data block
    % T = the total struct
    
    mask=T.dose.mask; % load the dose mask
    
    % pick the valid points from the mask 
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];
    
    fdose=zeros(1,size(d,3),'single'); % per-frame-dose
    fdosemin=zeros(1,size(d,3),'single'); % per-frame-minimal dose

    nframes=T.Raw.nframes(T.cas,T.rep); % parfor-hack
    

    parfor fra=1:nframes
        im=squeeze(d(:,:,fra));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fdose(fra)=mean(im)
        fdosemin(fra)=min(im(im>0));
    end
    
    % build a dose block and use it to normalize the data block
    doseblock=repmat(reshape(fdose,[1,1,length(fdose)]),[T.Raw.BS(1:2),1]);
    d_out=d./doseblock;
end

function DoseMaskCheck(im,T,fid)
    if T.proofs==1 % check the dose areas
        h=figure(fid);clf;
        set(gcf,'name','dose area mask check')
        
        subplot(1,3,1)
        imshow(im',[]);set(gca,'YDir','normal')
        title('original')
        
        subplot(1,3,2)
        imshow(T.dose.mask',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d mask',T.cas,T.rep))
        
        subplot(1,3,3)
        imshow(single(im.*(0.5+T.dose.mask/2))',[]);set(gca,'YDir','normal')
        title(f.FigTName('combined',fid,T.cas,T.rep))
        title('combined')
        fname=f.FigFileName('dose mask',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
        
    end
end

function DoseCorrectCheck1(T,dosemin,fid)
    if T.proofs==1
    % figure 5: dose & mindose check
        
        dosemean=squeeze(T.dose.doses(T.cas,T.rep,:));
        h=figure(fid);clf;

        subplot(2,1,1)
        plot(dosemean)
        title(f.FigTName('mean dose in dose areas per frame',fid,T.cas,T.rep))

        ylabel('mean dose')
        %xlabel('frame number')
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        grid on

        subplot(2,1,2)
        plot(dosemin)
        ylabel('minimum dose image value')
        xlabel('frame number')
        title(f.FigTName(sprintf(...
            'if this has major dips, then some object \n rotates into the dose mask'),...
            fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
    end
end

function DoseCorrectCheck2(d,T,fid)
    if T.proofs==1
    %check if the corrected values are ok
        h=figure(fid);clf;
        subplot(2,3,4:5);
        hold on
        ndf=size(T.dose.windows,2); % number of dose fields

        doseprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        %oldprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        dosepix=zeros(1,ndf);

        for i=1:ndf %iterate dose fields
            % mean of each field for each frame, for both old (im) and
            % corrected (imc) data
            doseprofile(:,i)=squeeze(mean(mean((d(...
                T.dose.windows(1,i):T.dose.windows(2,i),...
                T.dose.windows(3,i):T.dose.windows(4,i),:)))));
%             oldprofile(:,i)=squeeze(mean(mean((d(...
%                 T.dose.windows(1,i):T.dose.windows(2,i),...
%                 T.dose.windows(3,i):T.dose.windows(4,i),:)))));

            % plot
            plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))

            % the number of dose mask pixels
            dosepix(i)=(T.dose.windows(2,i)-T.dose.windows(1,i)+1)...
                *(T.dose.windows(4,i)-T.dose.windows(3,i)+1);

            % the total dose, should be constantly equal 1
            doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
%             oldprofile(:,ndf+1)=oldprofile(:,ndf+1)+dosepix(i)*oldprofile(:,i);
        end
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        %total dose manually done
        doseprofile(:,ndf+1)=doseprofile(:,ndf+1)/sum(dosepix);
%         oldprofile(:,ndf+1)=oldprofile(:,ndf+1)/sum(dosepix);
        plot(doseprofile(:,ndf+1),'DisplayName','combined')

        title(sprintf('relative doses of the %d regions \n after dose correction',ndf))
        xlabel('frame Nr')
        ylabel('relative dose')
        grid on

        legend()

        subplot(2,3,1:2);
        plot(squeeze(T.dose.doses(T.cas,T.rep,:)));
        ylabel('mean pixel value')
        title(f.FigTName('mean per frame dose of raw footage',fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])

        subplot(2,3,[3,6]);
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        hold on
        fcol=get(groot,'DefaultAxesColorOrder');
        for i=1:ndf
            r=rectangle('Position',[T.dose.windows(1,i),...
                T.dose.windows(3,i),...
                T.dose.windows(2,i)-T.dose.windows(1,i),...
                T.dose.windows(4,i)-T.dose.windows(3,i)],...
                'Facecolor',fcol(i,:));
        end
        title('dose areas')
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
    end
end


function [d_out,fits,fitdose,im_old]=BeamCorrect(d,T);
    % corrects for beam. now without dose correction
    % invariant for the dose according to the dose maks
    % d = data block
    % T = parameters struct
    % d_out = beam corrected data block
    % fits = cell array of the fit parameters
    % fitdose = post-beamcorrection-dose correction per frame
    % im_old is pre-beam-correction data for vizualisation purpose

    mask=T.dose.mask; %re-label
    
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single'); % contains the correction factor

    % find the fits
    %fprintf('%d %d finding beam non-uniformity fits...\n',T.cas,T.rep)
    nframes=T.Raw.nframes(T.cas,T.rep); % parfor-hack
    fitdose=zeros(1,nframes);
    im_old=squeeze(d(:,:,nframes));
    parfor j=1:nframes

        %f.f_BoCount(j,200,10,5)
        im=squeeze(d(:,:,j));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fits{j}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
        tbeam=feval(fits{j},x,y);         % make the beam
        fitdose(j)=mean(tbeam(mask==1));
        unibeam(:,:,j)=tbeam./fitdose(j);
   

    end
    % correct the non uniformity
    
    d_out=d./unibeam;

end

function BeamCorrectCheck(d,T,im_old)
    
    if T.proofs==1
        fid=60;
        h=figure(fid);clf;
        nframes=T.Raw.nframes(T.cas,T.rep);
        mask=T.dose.mask; %re-label
        
        % pick the valid points from the mask & prepare for fitting
        maskmesh=size(mask);
        [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
        xx=x;
        yy=y;
        xx(mask==0)=[]; % remove pixels outside of dose areas
        yy(mask==0)=[];
        im_old(mask==0)=[];
        imax=1;
        
        
        im=squeeze(d(:,:,nframes));
        im(T.dose.mask==0)=[];
        rf =10; % reduction factor
       
        
        plot(T.fit.fits{T.cas,T.rep}{nframes},[xx(1:rf:end)',yy(1:rf:end)'],...
            im_old(1:rf:end)')
        hold on
        tfit= fit( [xx',yy'],double(im'), 'poly22' );
        tbeam=feval(tfit,x,y); 
        surf(x(1:rf:end,1:rf:end),y(1:rf:end,1:rf:end),...
            tbeam(1:rf:end,1:rf:end))
        
        %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
        xlabel('x')
        ylabel('y')
        zlabel('counts')
        %pause(0.5)
        
        view([20,45])
        title(f.FigTName(sprintf(...
            'beam nun-uniformity fit \n and corrected image'),...
            fid,T.cas,T.rep))
        fname=f.FigFileName('Beam fit',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
    end
end

function do=dose(im,mask)
    % returns the dose according to the mask of an image
    im(mask==0)=[];
    do=mean(im);
end

function [fitim,fits]=FrameFit(im,mask)
    % returns the fit of frame according to a mask
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];
    im(mask==0)=[];
    fits= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
    fitim=feval(fits,x,y);
end

function [d_out,fi,unibeam]=BeamAndDoseCorrection3(d,T)
    % corrects for beam & dose
    % d = data block
    % T = parameters struct
    %fprintf('%d %d correctiong for beam and dose...\n',T.cas,T.rep)
    mask=T.dose.mask; %re-label
    
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single'); % contains the correction factor
    fdose=zeros(size(d,3)); % per-frame-dose

    % find the fits
    %fprintf('%d %d finding beam non-uniformity fits...\n',T.cas,T.rep)
    nframes=T.Raw.nframes(T.cas,T.rep); % parfor-hack
    parfor j=1:nframes

        %f.f_BoCount(j,200,10,5)
        im=squeeze(d(:,:,j));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fdose(j)=mean(im)
        fi{j}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
        unibeam(:,:,j)=feval(fi{j},x,y);         % make the beam

    end
   
    %fprintf('\n')

    % correct the non uniformity
    d_out=d./unibeam;

    % Dose Area checks
    if T.proofs==1 % check the areas
        fid=40;
        h=figure(fid);clf;
        set(gcf,'name','dose area mask check')
        
        subplot(1,3,1)
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        title('original')
        
        subplot(1,3,2)
        imshow(T.dose.mask',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d mask',T.cas,T.rep))
        
        subplot(1,3,3)
        imshow((single(squeeze(d(:,:,1))).*(0.5+T.dose.mask/2))',[]);set(gca,'YDir','normal')
        title(f.FigTName('combined',fid,T.cas,T.rep))
        title('combined')
        fname=f.FigFileName('dose mask',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)



        %check if nothing rotates in these areas
        fid=50;
        h=figure(fid);clf;
        dosemin=zeros(1,T.Raw.nframes(T.cas,T.rep));
        dosemean=zeros(1,T.Raw.nframes(T.cas,T.rep));
        for i=1:T.Raw.nframes(T.cas,T.rep)
            a=single(squeeze(d(:,:,i))).*T.dose.mask;
            dosemin(i)=min(a(a>0));
            dosemean(i)=mean(a(a>0));
        end

        subplot(2,1,1)
        plot(dosemean)
        title(f.FigTName('mean dose in dose areas per frame',fid,T.cas,T.rep))

        ylabel('mean dose')
        xlabel('frame number')
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        grid on

        subplot(2,1,2)
        plot(dosemin)
        ylabel('minimum dose image value')
        xlabel('frame number')
        title(f.FigTName(sprintf(...
            'if this has major dips, then some object \n rotates into the dose mask'),...
            fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)


    end

    %fit checks
    if T.proofs==1
        fid=60;
        h=figure(fid);clf;
        imax=1;
        ii=ceil(T.Raw.nframes(T.cas,T.rep)*rand(1,imax)); %check some random frames
        for i = 1:length(ii)

            im=squeeze(d(:,:,ii(i)));
            im(T.dose.mask==0)=[];
            rf =10; % reduction factor
            plot(fi{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
            hold on
            %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
            xlabel('x')
            ylabel('y')
            zlabel('counts')
            %pause(0.5)
        end
        view([20,45])
        title(f.FigTName(sprintf(...
            'nun-uniformity data \n according to a random frame'),...
            fid,T.cas,T.rep))
        fname=f.FigFileName('Beam fit',fid,T.cas,T.rep);    
        f.f_BoFig2PDF(h,fname)

        %check if the corrected values are ok
        fid=70;
        h=figure(fid);clf;
        subplot(2,3,4:5);
        hold on
        ndf=size(T.dose.windows,2); % number of dose fields

        doseprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        oldprofile=zeros(T.Raw.nframes(T.cas,T.rep),ndf+1);
        dosepix=zeros(1,ndf);

        for i=1:ndf %iterate dose fields
            % mean of each field for each frame, for both old (im) and
            % corrected (imc) data
            doseprofile(:,i)=squeeze(mean(mean((d_out(...
                T.dose.windows(1,i):T.dose.windows(2,i),...
                T.dose.windows(3,i):T.dose.windows(4,i),:)))));
            oldprofile(:,i)=squeeze(mean(mean((d(...
                T.dose.windows(1,i):T.dose.windows(2,i),...
                T.dose.windows(3,i):T.dose.windows(4,i),:)))));

            % plot
            plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))

            % the number of dose mask pixels
            dosepix(i)=(T.dose.windows(2,i)-T.dose.windows(1,i)+1)...
                *(T.dose.windows(4,i)-T.dose.windows(3,i)+1);

            % the total dose, should be constantly equal 1
            doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
            oldprofile(:,ndf+1)=oldprofile(:,ndf+1)+dosepix(i)*oldprofile(:,i);
        end
        xlim([1,T.Raw.nframes(T.cas,T.rep)])
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
        title(f.FigTName('mean per frame dose of raw footage',fid,T.cas,T.rep))
        grid on
        xlim([1,T.Raw.nframes(T.cas,T.rep)])

        subplot(2,3,[3,6]);
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        hold on
        for i=1:ndf
            r=rectangle('Position',[T.dose.windows(1,i),...
                T.dose.windows(3,i),...
                T.dose.windows(2,i)-T.dose.windows(1,i),...
                T.dose.windows(4,i)-T.dose.windows(3,i)],...
                'Facecolor','r');
        end
        title('dose areas')
        fname=f.FigFileName('dose area checks',fid,T.cas,T.rep);
        f.f_BoFig2PDF(h,fname)
        


    end
   
end

function [d_out]=SimpleBeamCorrection(d,T)
   % simply applies the masterbeam to the data
    % d = data block
    % T = parameters struct

    mask=T.dose.mask; %re-label
    
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single'); % contains the correction factor

    % find the fits
    %fprintf('%d %d finding beam non-uniformity fits...\n',T.cas,T.rep)
    nframes=T.Raw.nframes(T.cas,T.rep); % parfor-hack
    masterbeam=T.fit.masterbeam;
    for j=1:nframes
        im=squeeze(d(:,:,j));    % check out a slice from stack
       meandose=mean(im(mask));
       unibeam(:,:,j)=masterbeam.*meandose;         % make the beam
    end
   
   d_out=d./unibeam;
end

function T=PreProcessingPreAllocation(T)
    % returns empty arrays. 
    % main pupose is to keep the main script tidy
    fprintf('Pre-prcessing...')
    T.Raw.rad=zeros([T.d.imsize,T.d.ncas,T.d.nrep],T.sys.StoreFormat);    % radiography
    T.Raw.CropRange=zeros(T.d.ncas,T.d.nrep,2,3); % stores the crop range data
    T.Raw.nframes=zeros(T.d.ncas,T.d.nrep); %nframes after cropping
    T.Match.shift=zeros(T.d.ncas,T.d.nrep,3); % XYZ shifts from matching
    T.Match.diff=zeros(T.d.ncas,T.d.nrep,3);  % XZ shift differences
    T.dose.mask=f.MakeMask(T.dose.windows,T.d.imsize); %creates the dosemask
    
    
    % load the first sinogram for once before the big loop to simplify
    % stuff
    cas=1;
    rep=1;
    T.cas=1;
    T.rep=1;
    
    % 1 read in
    d=f.ReadXrayTomo(T);
    
    % 2 find rot startt & stop and crop
            croprange=f.FindCrop(T,...
                squeeze(d(:,T.Raw.sinoyplanes(1),:)));
    T.Raw.CropRange(cas,rep,:,1)=croprange(:,1);
    T.Raw.CropRange(cas,rep,:,2)=croprange(:,1);
    [d,BS3]=f.SinoCrop(d,T);
    BS=[T.d.imsize(1),T.d.imsize(2),BS3];
    T.Raw.BS=BS;
    T.Raw.Sino=zeros(BS(1),BS(3),T.d.ncas,T.d.nrep,size(T.Raw.sinoyplanes,2),'single');
           
    %creates the Shiftmask xy
    T.Match.MaskSino=f.MakeMask(T.Match.SinoWindows,...
        [BS(1),BS(3)]);
    T.Match.MaskY=f.MakeMask(T.Match.YWindows,...
        [BS(2),BS(3)]);
      
    % q360 data
    T.q360.nFrames=zeros(T.d.ncas,T.d.nrep);
    T.q360.nFrames(T.cas,T.rep)=f.frames360(squeeze(d(:,T.Raw.sinoyplanes(1),:)),T);
    T.Log.Sino=zeros(T.Raw.BS(1),T.q360.nFrames(1,1),T.d.ncas,'single'); % divved sinos
    T.fnames.add=cell(1,T.d.ncas);
    T.fnames.dive=cell(1,T.d.ncas);
    T.fnames.corr=cell(T.d.ncas,T.d.nrep);
    T.fnames.corr2=cell(T.d.ncas,T.d.nrep);
    T.fnames.rec=cell(T.d.ncas,T.d.nrep);
    
    fprintf('done. \n')
    
    % load the masterbeam, the "mean shape" of the beam
    T.fit.masterbeam=f.BoLoad('masterbeam',T);
    
    % dose data block
    T.dose.doses=zeros(T.d.ncas,T.d.nrep,BS3,'single');
    
    T.Cen.eshift=zeros(T.Raw.BS(2),T.d.ncas,T.d.nrep,'single');
    T.Cen.xlist=zeros(T.Cen.niter,T.Raw.BS(2),T.d.ncas,T.d.nrep,'single');
    T.Cen.v=nan(T.Cen.niter,T.Raw.BS(2),T.d.ncas,T.d.nrep,'single');
    T.Cen.MMask=zeros(length(T.Cen.range),length(T.Cen.range)...
        ,T.d.ncas,T.d.nrep,'single'); 
    T.Fig=f.FigProperties(T);
end

function CheckSinoShiftMask(T)
    %% manually find the shift diff area
    % just plots the mask
    h=figure(106);clf
    im=squeeze(T.Raw.Sino(:,:,1,1,1)).*(0.5*T.Match.DiffMaskSino+0.5);
    imshow(im,[])
    fname=sprintf('shift base area plot');
    f.f_BoFig2PDF(h,fname)
end

function s=FigFileName(text,fid,cas,rep)
    s=sprintf('Fig %03d %02d %02d %s',fid,cas,rep,text);
end

function s=FigTName(text,fid,cas,rep)
    s=sprintf('F%d [%d,%d]: %s',fid,cas,rep,text);
end

function s=FigFileName2(text,fid,cas)
    s=sprintf('Fig %03d %02d %s',fid,cas,text);
end

function s=FigTName2(text,fid,cas)
    s=sprintf('F%d [%d]: %s',fid,cas,text);
end

function s=FigFileName3(text,fid,cas,rep,plane)
    s=sprintf('Fig %03d %02d %02d %04d %s',fid,cas,rep,plane,text);
end



%fname=f.FigFileName2(text,fid,T.cas)
%title(f.FigTName2(text,fid,cas))
function s2=CommLineStart(cas,rep)
    % prepares a unified string for command lineoutputs
    % 14:34 5h 3min 05 09 
    s=datestr(now(),'hh:MM');
    try
        b=toc;
    catch
        tic;
        b=toc;
    end
        if b<60
       c= sprintf('%2.0fs',b);
    elseif b<3600
       c= sprintf('%2.0fm',b/60);
    else
         c= sprintf('%2.0fh',b/3600);
    end
    s2=sprintf('%s (+%s) [%02d %02d]',s,c,cas,rep);
end

function s2=CommLineStart2(cas)
    % prepares a unified string for command lineoutputs
    % 14:34 5h 3min 05 09 
    s=datestr(now(),'hh:MM');
    try
        b=toc;
    catch
        tic;
        b=toc;
    end
        if b<60
       c= sprintf('%2.0fs',b);
    elseif b<3600
       c= sprintf('%2.0fm',b/60);
    else
         c= sprintf('%2.1fh',b/3600);
    end
    s2=sprintf('%s (+%s) [%02d]',s,c,cas);
end

function corr=SinglePixCorr(im,sp)
    % corrects the few pixel errors
    % base image: mean of all projectuions
    corr=im;
    
    % 1) single pixels. idea: replace faulty pixels by the mean of the 4
    % pixels to the north, south, east, west
    
    modx=[1,-1,0,0];
    mody=[0,0,1,-1];
    
    idx=repmat(sp(:,1),[1,4])+repmat(modx,[size(sp,1),1]);
    idy=repmat(sp(:,2),[1,4])+repmat(mody,[size(sp,1),1]);
    
    corrpix=mean(im( sub2ind(size(im), idx, idy)),2);
    corr(sub2ind(size(corr),sp(:,1),sp(:,2)))=corrpix;
    
end    
    
function corr=SpecialPixFix(im)
    
    corr=im;
    % special
    corr(262,355)=mean([corr(262,355+1),corr(262,355-1),corr(262-1,355)]);
    
    % linear faults corected via x-mean
    x=262;
    y=79:81;
    corr(x,y)=mean([corr(x-1,y);corr(x+1,y)],1);
    
    x=263;
    y=349:360;
    corr(x,y)=mean([corr(x-1,y);corr(x+1,y)],1);
    
    x=333;
    y=854:855;
    corr(x,y)=mean([corr(x-1,y);corr(x+1,y)],1);
    
    
    % linear faults corected via y-mean
    x=409:411;
    y=588;
    corr(x,y)=mean([corr(x,y-1);corr(x,y+1)],1);
    
    x=309:311;
    y=139;
    corr(x,y)=mean([corr(x,y-1);corr(x,y+1)],1);
    
end

function corr_out=ImCorrFilt(im,T)
    % fixes some camera spots and applies a soft gauss filter
    % im: mean(d,3)
    % T: T struct
    % corr: corrected image
    
    % correct the zero values, usually 5 rows typically at start & end (because
    % they fuck with the conv2 filter later)
    zthresh=200; % from befor d normalization
    zthresh=0.5;
    z=(mean(im)<zthresh);
    last=find(z(1:round(T.Raw.BS(2)/2))==1,1,'last');
    z(1:last)=0;
    first=find(z==1,1,'first');
    im(:,1:last+1)=repmat(im(:,last+2),[1,last+1]);
    im(:,first-1:end)=repmat(im(:,first-2),[1,T.Raw.BS(2)-first+2]);
    % and fix first row
    im(1,1:136)=im(2,1:136);

    
    % 1) camera spot removal: singel pixels
    sp=[264, 468;... % faulty single pixels x,y;
        254, 388;... single pixel
        411, 510;... single pixel
        315, 752;... single pixel
        252, 760;... single pixel
        239, 744;... single pixel
        262, 718;... single pixel
        270, 921;... single pixel
        382, 770;...
        332, 781;...
        327, 410;...
        342, 895;...
        276, 341;...
        414, 525;...
        400, 538;...
        329, 611;...
        390, 533;...
        259, 780;...
        236, 72;...
        312, 464;
        385, 377;
        333, 685;
        326, 726;
        333, 935;
        354, 964;
        261, 926;
        206, 681;
        294, 473];
    
    %apply this list:
    im2=f.SinglePixCorr(im,sp);
    
    % special faults:
    im3=f.SpecialPixFix(im2);
    
    % blob faults:
    im4=f.RemoveCameraSpots(im3,T.cas,T.rep);
    
    % that weired structure that pissed me off:
    xw1=322:333;
    yw1=573*ones(1,length(xw1));
    xw2=319:328;
    yw2=574*ones(1,length(xw2));
    xw3=313:322;
    yw3=575*ones(1,length(xw3));
    xw4=300:305;
    yw4=579*ones(1,length(xw4));    
    xw5=302:306;
    yw5=578*ones(1,length(xw5));
    xw6=303:311;
    yw6=577*ones(1,length(xw6));
    xw7=310:314;
    yw7=576*ones(1,length(xw7));
    xw8=326:333;
    yw8=572*ones(1,length(xw8));
    xw9=[297,297];
    yw9=[580,581];
    xw10=[293:294];
    yw10=[582,582];
    xw11=290:294;
    yw11=583*ones(1,length(xw11));
    xw12=[291,284,282,278,277,277,274,272,270,265];
    yw12=[584,587,589,591,591,592,593,594,596,599];
    
    wx=[xw1,xw2,xw3,xw4,xw5,xw6,xw7,xw8,xw9,xw10,xw11,xw12];  % "mask" for that structure
    wy=[yw1,yw2,yw3,yw4,yw5,yw6,yw7,yw8,yw9,yw10,yw11,yw12];
    
    im5=f.RemoveAnnoyingStructre(im4,wx,wy);
    
    corr=im5;
    % apply gauss filter. mor in y- and less in x direction
    % convolute an elliptical gaussian Kernel with the image
    % i try a few and later chose one
    
    kern=mt.ElliKernelGen(3,7,0.5,1);
    kern=mt.ElliKernelGen(1,7,1,1);
    im6=imfilter(im5,kern,'replicate');

    kern=mt.ElliKernelGen(3,9,0.5,2);
    im7=imfilter(im5,kern,'replicate');
    
    kern=mt.ElliKernelGen(3,13,0.5,3);
    im8=imfilter(im5,kern,'replicate');
    
    kern=mt.ElliKernelGen(3,17,0.5,4);
    im9=imfilter(im5,kern,'replicate');    
    
    kern=mt.ElliKernelGen(3,25,0.5,5);
    im10=imfilter(im5,kern,'replicate');  
    

    corr(:,401:550)=im6(:,401:550);
    corr(:,551:700)=im7(:,551:700);
    corr(:,701:950)=im8(:,701:950);
    corr(:,701:950)=im9(:,701:950);
    corr(:,951:end)=im9(:,951:end);
    
    %decision: take the mild gauss filter
    corr=im6;
    
    corr_out=corr./im;
    
    crange=[2.5,3.1]*10000; % from befor normalizing d
    crange=[0.55,0.7]; % from befor normalizing d
    
    imrange=200:450;
    
    fid=90;
    h=figure(fid);clf;
    nplots=5;
    subplot(1,nplots,1)
    imagesc(im(imrange,:)');
    title(f.FigTName('original',fid,T.cas,T.rep));
    ax(1)=gca;
    colormap(hsv)
    caxis(crange);
    set(gca,'YDir','normal')
    axis off
    axis equal
    axis tight
    
    subplot(1,nplots,2)
    imshow((im5(imrange,:)~=im(imrange,:))');
    title('fixed pixels')
    ax(2)=gca;
    set(gca,'YDir','normal')
    axis off
    axis equal
    axis tight
    
    subplot(1,nplots,3)
    imagesc(im5(imrange,:)');
    title('pixels corrected')
    ax(3)=gca;
    colormap(hsv)
    caxis(crange);
    set(gca,'YDir','normal')
    axis off
    axis equal
    axis tight
    
    subplot(1,nplots,4)
    imagesc(corr(imrange,:)');
    title('gauss filtered')
    ax(4)=gca;
    colormap(hsv)
    caxis(crange);
    set(gca,'YDir','normal')
    axis off
    axis equal
    axis tight
        
    subplot(1,nplots,5)
    imshow(corr_out(imrange,:)');
    caxis(1+[-0.01,+0.01])
    title('relative fix')
    ax(5)=gca;
    set(gca,'YDir','normal')
    axis off
    axis equal
    axis tight
    
    linkaxes(ax)
    fname=f.FigFileName('Image filter',fid,T.cas,T.rep);
    f.f_BoFig2PDF(h,fname)
    
end

function corr=RemoveAnnoyingStructre(im,wx,wy)
    % this function removes the detector artefact by overwriting the pixels
    % specified in [wx,wy] by an average value of 3 pixels above and below
    % (the structure is somewhat horizontal)
    corr=im;
    ux=unique(wx);
    for x=ux
        ymax=max(wy((wx-x)==0));
        ymin=min(wy((wx-x)==0));
        corr(x,ymin:ymax)=mean([im(x,(ymin-3):(ymin-1)),im(x,(ymax+1):(ymax+3))]);
    end
end

function d=ApplyImCorr(d,corr,T)
    % applyies the pixel fixes and filter to the d block
    d2=repmat(corr,[1,1,size(d,3)]);
    d=d.*d2;
end

function structstruct(S)

% Figure the type and class of the input
whosout = whos('S');
sizes = whosout.size;
sizestr = [int2str(sizes(1)),'x',int2str(sizes(2))];
by=whosout.bytes;
if by <1024
    bystr=sprintf(' %d B',by);
elseif by<1024*1024
    bystr=sprintf(' %.2f kB',by/1024);
elseif by<1024*1024*1024
    bystr=sprintf(' %.2f MB',by/1024/1024);
else
    bystr=sprintf(' %.2f GB',by/1024/1024/1024);
end

endstr = [':  [' sizestr '] ' whosout.class bystr];

% Print out the properties of the input variable
disp(' ');
disp([inputname(1) endstr]);

% Check if S is a structure, then call the recursive function
if isstruct(S)
    f.recursor(S,0,'');
end

% Print out a blank line
disp(' ');

end

function recursor(S,level,recstr)

recstr = [recstr '  |'];

fnames = fieldnames(S);

for i = 1:length(fnames)
    
    %% Print out the current fieldname
    
    % Take out the i'th field
    tmpstruct = S.(fnames{i});
    
    % Figure the type and class of the current field
    whosout = whos('tmpstruct');
    sizes = whosout.size;
    ss=size(sizes);
    sizestr = int2str(sizes(1));
    for j=2:length(sizes);
        sizestr=strcat(sizestr,'x',num2str(sizes(j))); 
    end
    sizestr=['[' sizestr ']'];
    
    by=whosout.bytes;
    if by <1024
        bystr=sprintf(' %3.0f     B',by);
    elseif by<1024*1024
        bystr=sprintf(' %3.2f kB',by/1024);
    elseif by<1024*1024*1024
        bystr=sprintf(' %3.2f MB',by/1024/1024);
    else
        bystr=sprintf(' %3.2f GB',by/1024/1024/1024);
    end
    
    endstr = [':  [' sizestr '] ' whosout.class bystr];
    endstr =sprintf(': %17s %6s %10s ',sizestr,whosout.class, bystr);
    
    % Create the strings
    if i == length(fnames) % Last field in the current level
        startstr=[recstr(1:(end-1)) '''--' fnames{i}];
       recstr(end) = ' ';
    else % Not the last field in the current level
        startstr=[recstr '--' fnames{i}];
    end
    str=sprintf('%-20s %s',startstr, endstr);
    % Print the output string to the command line
    disp(str);
    
    %% Determine if each field is a struct
    
    % Check if the i'th field of S is a struct
    if isstruct(tmpstruct) % If tmpstruct is a struct, recursive function call
        f.recursor(tmpstruct,level+1,recstr); % Call self
    end
    
end

end

function GetSize(this)
    props = properties(this);
    totSize = 0;
    for ii=1:length(props)
        currentProperty = getfield(this, char(props(ii)));
        s = whos('currentProperty');
        totSize = totSize + s.bytes;
    end
    fprintf(1, '%d bytes\n', totSize);
end

function F=FigProperties(T)
    % defines centrally figure properties
% x labels
F.xl.FU='points';   % Font size unit
F.xl.FS=11;         % Font size
F.xl.FN='Times';    % Font name
F.xl.Interp='latex';% Latex interpreter
F.xl.FW='normal';    % Font weight

% y labels
F.yl.FU='points';   % Font size unit
F.yl.FS=11;         % Font size
F.yl.FN='Times';    % Font name
F.yl.Interp='latex';% Latex interpreter
F.yl.FW='normal';    % Font weight

% title
F.ti.FU='points';   % Font size unit
F.ti.FS=12;         % Font size
F.ti.FN='Times';    % Font name
F.ti.Interp='latex';% Latex interpreter
F.ti.FW='bold';      % Font weight

% axes
F.ax.FU='points';   % Font size unit
F.ax.FS=10;         % Font size
F.ax.FN='Times';     % Font name
F.ax.Interp='latex';% Latex interpreter
F.ax.FW='normal';    % Font weight
F.ax.U='normalized'; % ???

% figure properties 
F.Fig.Units='centimeters';  
F.Fig.PPM='auto';       %PaperPositionMode    
F.Fig.FigW=8; % determines figure sizes
F.Fig.FigH=6; % both taken form the latex textwidth & 4:3 aspec tatio
F.Fig.PapS=[8,6];

%Legend
F.L.FU='points';   % Font size unit
F.L.FS=8;         % Font size
F.L.FN='Times';     % Font name
F.L.Interp='latex';% Latex interpreter
F.L.FW='normal';    % Font weight



% line plots
F.L.LW=2;            %LineWidth

F.saveto='V:\1phd\thesis\figures\'; % save folder

end

function [shift,xlist,v]=fibo(sino,baseshift,xmin,xmax,imax,epsilon,key,ang,...
        recsize,detPitch,src,det,mask)
    % performs the fibonacci center finding
    tau=double((sqrt(5)-1)/2);
    x1n=@(aa,b,tau)aa+(1-tau)*(b-aa);   % lefz mid point
    x2n=@(aa,b,tau)aa+tau*(b-aa);      % riht mid point
    
    %preallocate
    xlist=nan(1,imax);
    v=nan(1,imax);
    rec=zeros(recsize,recsize,imax);
    
    % fibonacci spacing
    aa=xmin; % intervall bounds
    b=xmax;
    i=2;  % number of (initial) iterations, counter
    
    x1=x1n(aa,b,tau);             % computing initial mid x values
    x2=x2n(aa,b,tau);
    xlist(1:2)=[x1,x2]; % list for plotting
    
    % find the two initial values

    if key==1 % parallel
        [f1,~]=f.centerReconPar(sino,baseshift+x1,ang,recsize,mask);
        [f2,~]=f.centerReconPar(sino,baseshift+x2,ang,recsize,mask);
    
    elseif key==2 % fan
        [f1,~]=f.centerReconFan(sino,baseshift+x1,ang,recsize,...
            detPitch,src,det,mask);
        [f2,~]=f.centerReconFan(sino,baseshift+x2,ang,recsize,...
            detPitch,src,det,mask);
    end    
    
    v(1:2)=[f1,f2];

    while ((abs(b-aa)>epsilon) && (i<=imax))
        i=i+1;
        if(f1>f2)
            b=x2;
            x2=x1;
            x1=x1n(aa,b,tau);
            xlist(i)=x1;
            f2=f1;
            
            if key==1 % parallel
                [f1,~]=f.centerReconPar(sino,baseshift+x1,ang,recsize,mask);
                
            elseif key==2 % fan
                [f1,~]=f.centerReconFan(sino,baseshift+x1,ang,recsize,...
                    detPitch,src,det,mask);
            end
            v(i)=f1;
        else
            aa=x1;
            x1=x2;
            x2=x2n(aa,b,tau);
            xlist(i)=x2;
            
            f1=f2;
            
            if key==1 % parallel
                [f2,~]=f.centerReconPar(sino,baseshift+x2,ang,recsize,mask);
                
            elseif key==2 % fan
                [f2,~]=f.centerReconFan(sino,baseshift+x2,ang,recsize,...
                    detPitch,src,det,mask);
            end
            v(i)=f2;
        end
    end
    
    [vmax,vind]=max(v);
    shift=xlist(vind);
end

function [q,rec]=centerReconFan(sino,shift,ang,recsize,...
        detPitch,src,det,FanMask)
    sinogram=f.fraccircshift(sino,...
        -shift);
    rec=a.FBPexplFan(sinogram',...
        recsize,ang,detPitch,src,det);
    [q,~]=qc.VarianceQuality(rec,FanMask);
end

function [q,rec]=centerReconPar(sino,shift,ang,recsize,ParMask)
    sinogram=f.fraccircshift(sino,...
        -shift);
    rec=a.FBPexpl(sinogram',recsize,ang);
    [q,~]=qc.VarianceQuality(rec,ParMask);
end

function mask=MakeTomoMask(key,sino,recsize,ang,detPitch,src,det,...
        threshpar,threshfan,point)
    % make a mask for the edge matching.
    % sino =  sinogram
    % recsize
    % ang
    % thresh for the imbw operation
    % point for filling ( remove film edge from mask)
    
    % make first mask
    
    if ischar(key)
        if key =='par'
            recmask=a.FBPexpl(sino',recsize,ang);
            [~,mask,~]=qc.kern(recmask,threshpar,point);
        elseif key == 'fan'
            recmask=a.FBPexplFan(sino',recsize,ang,detPitch,src,det);
            [~,mask,~]=qc.kern(recmask,threshfan,point);
        else
            error('mask key must be "par" or "fan"')
        end
    else
        error('key input variable must be string: "par" or "fan"')
    end
        
end



function CheckRecon(rec,v,xlist,cas,rep,plane,xmin,xmax,writekey)
    fid=140;
    fig=figure(fid);clf;
    fig.Position=[100 162 400 800];
    subplot('Position',[0.05,0.5,0.9,0.45]);
    imshow(rec,[]);
    title(sprintf('c%d r%d p%03d',cas,rep,plane));
    subplot('Position',[0.1,0.1,0.8,0.3]);
    
    pv=v;
    pind=~isnan(pv);
    plot(xlist(pind),v(pind),'xr','Displayname','fibo-algo')
    hold on
    [vmax,vind]=max(v);
    
    plot(xlist(vind),v(vind),'+b','Displayname',...
        sprintf('%.3f',xlist(vind)))
    title(f.FigTName('Another mask',fid,cas,rep))
    title(sprintf('shift finding, %d iterations',length(pind)))
    grid on
    xlim([xmin,xmax]);
    legend('Location','South')
    
    if writekey
        fname=f.FigFileName3('ReconCheck',fid,cas,rep,plane);
        f.f_BoFig3PDF(fig,'reconcheck\',fname)
    end
    
end
    
function CheckMMask(MMask,cas,rep)
    fid=135;
    fig=figure(fid);clf;
    imshow(MMask,[]);
    title(f.FigTName('Another mask',fid,cas,rep))
    fname=f.FigFileName('MMaskCheck',fid,cas,rep);
    f.f_BoFig2PDF(fig,fname)

end

function [list,cas,plane,raw]=GenRecFileList(T,fid)
    % read the file names and preallocate the data
    oldPath=cd(T.pp.loadpath);
    list=dir('Micha*.csv');
    cd(oldPath);
    cas=zeros(1,length(list));
    plane=zeros(1,length(list));
    
    % read in every file, extract case & plane
    for i=1:length(list)
        name=list(i).name;
        cas(i)=str2num(name(T.pp.casstr));
        plane(i)=str2num(name(T.pp.plastr));
    end
    
    % sum up
    raw.ncas=length(unique(cas));
    raw.cases=unique(cas);
    raw.planes=unique(plane);   
    
    % check how many of each
    for i = 1:raw.ncas
        raw.nplane(i)=sum(cas==raw.cases(i));
        
    end
    
    
    
end

function mtomo=ReadMichasTomos(F,rs,fid)
    % reads in the micha tomos
    oldPath=cd(F.pp.loadpath);
    mtomo=nan(rs,rs,F.raw.ncas,length(F.planes),'single');
    fprintf('reading in mtomo....')
    
    for i=1:length(F.pp.flist)
        name=F.pp.flist(i).name;
        cas=str2num(name(F.pp.casstr));
        plane=str2num(name(F.pp.plastr));
        plid=find(F.planes==plane);
        casid=find(F.raw.cases==cas);
        mtomo(:,:,casid,plid)=dlmread(F.pp.flist(i).name);
        
    end
    cd(oldPath);
    fprintf('done. saving mtomo.....')
    save('mtomo','mtomo')
    fprintf('done.\n')
    
end

function planes=workplaneslist(F)
    % since michas algo does not converge for the spacer planes, we have to
    % artificially insert them here for consisten plotting
    
    
    
end

%
function sumrange=LFTSum(F,T) % find the best sumrange
end

function im=imnorm(im)
% renormalize an image
im=(im-min(im(:)))/(max(im(:))-min(im(:)));
end


end %static
end %class





















