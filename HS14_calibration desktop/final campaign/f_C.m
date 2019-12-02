classdef f_C < handle


properties
    T
end
%% function collection file for the final campaign Dec. 201+
methods(Static) %evil function scope hack
    function obj=f_C()
        obj.T = []
    end
    
    function out=ClassTest()
        disp('hello world')
        out.T
    end

function obj=GenerateFileList(obj)
    % generates the files list in the folder, along with some parameters
    reps=zeros(obj.T.d.ncas,obj.T.d.nrep); % indicator of available files
    oldPath=cd(obj.T.d.DataPath);
    list=dir('*tom_*.seq');
    cd(oldPath);
    filenames=cell(obj.T.d.ncas,obj.T.d.nrep);
    nframes=zeros(obj.T.d.ncas,obj.T.d.nrep);
    for i=1:length(list) % iterate mesurements
        fname=list(i).name;
        cas=str2num(fname(14:15));
        rep=str2num(fname(17:18));
        filenames{cas,rep}=list(i).name;
        nframes(cas,rep)=(list(i).bytes-obj.T.d.header)...
            /obj.T.d.GreyBytes/prod(obj.T.d.imsize);
        if mod(nframes(cas,rep),1) ~= 0
            error('file size is not coherent with pixel sizes and header');
        end
        obj.T.d.reps(cas,rep,1)=1;
    end
    obj.T.d.nframes=nframes;
    obj.T.d.list=list;
    
    
    fid=13;
    h=figure(fid);clf
    d=nframes;
    d(d~=0)=d(d~=0)-min(d(d~=0));
    imagesc(d');colorbar()
    ax=gca;
    title('difference of tomo lenghts (in frames)')
    ylabel('repetition')
    xlabel('case/ betriebspunkt')
    zlabel('additional frames')
    fname=sprintf('Fig %02d Lengths differeces of raw tomos' ,fid);
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
    
    h.Renderer='Painters';
    
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

    fpath=strcat(obj.T.d.DataPath,obj.T.d.List{obj.T.cas,obj.T.rep});
    fprintf('reading in %s... \n',obj.T.d.List{obj.T.cas,obj.T.rep})
    % open the file
    fid=fopen(fpath);

    % dispose of header
    fread(fid,obj.T.d.header);%the header size migth need to be adjusted depending on image settings

    %read in data
    d=(fread(fid,obj.T.d.imsize(1)*obj.T.d.imsize(2)*obj.T.d.nframes(obj.T.cas,obj.T.rep),'uint16=>single')); %the image
    fclose(fid); %close file handle

    % rotate such that the representation is correct
    d=permute(reshape(d,[obj.T.d.imsize(2) obj.T.d.imsize(1) obj.T.d.nframes(obj.T.cas,obj.T.rep)]),[2,1,3]);

end

function Dosemask=MakeDoseMask(T)
    % creates the dose mask
    Dosemask=zeros(obj.T.d.imsize); %make the dose mask
    for i=1:size(obj.T.dose.windows,2)
        Dosemask(obj.T.dose.windows(1,i):obj.T.dose.windows(2,i),...
               obj.T.dose.windows(3,i):obj.T.dose.windows(4,i))=1;
    end
end

function mask=MakeMask(windows,imsize)
    % creates the dose mask
    mask=zeros(imsize); %make the dose mask
    for i=1:size(windows,2)
        mask(windows(1,i):windows(2,i),...
               windows(3,i):windows(4,i))=1;
    end
end

function [d_out,T]=BeamAndDoseCorrection(d,T)
    % corrects for beam & dose
    % d = data block
    % T = parameters struct

    mask=obj.T.mask; %re-label
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single');

    % find the fits
    fprintf('%d finding beam non-uniformity fits...\n',obj.T.i)
    nframes=obj.T.d.nframes(obj.T.cas,obj.T.rep); % parfor-hack
    parfor j=1:nframes

        f.f_BoCount(j,50,10,5)
        im=squeeze(d(:,:,j));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fi{j}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
        unibeam(:,:,j)=feval(fi{j},x,y);         % make the beam

    end
    obj.T.f=fi;
    fprintf('\n')

    % correct the non uniformity
    d_out=d./unibeam;
    obj.T.unibeam=unibeam;

    % Dose Area checks
    if obj.T.proofs==1 % check the areas

        h=figure(4);clf;
        set(gcf,'name','dose area mask check')
        subplot(1,3,1)
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        title(sprintf('%02d original',obj.T.i))
        subplot(1,3,2)
        imshow(obj.T.mask',[]);set(gca,'YDir','normal')
        title(sprintf('%02d mask',obj.T.i))
        subplot(1,3,3)
        imshow((single(squeeze(d(:,:,1))).*(0.5+obj.T.mask/2))',[]);set(gca,'YDir','normal')
        title(sprintf('%02d combined',obj.T.i))
        fname=sprintf('Fig04 %02d dose mask',obj.T.i);
        f.f_BoFig2PDF(h,fname)

        %check if nothing rotates in these areas
        h=figure(5);clf;
        dosemin=zeros(1,obj.T.d.nframes(obj.T.cas,obj.T.rep));
        dosemean=zeros(1,obj.T.d.nframes(obj.T.cas,obj.T.rep));
        for i=1:obj.T.d.nframes(obj.T.cas,obj.T.rep)
            a=single(squeeze(d(:,:,i))).*obj.T.mask;
            dosemin(i)=min(a(a>0));
            dosemean(i)=mean(a(a>0));
        end

        subplot(2,1,1)
        plot(dosemean)
        title(sprintf('%02d mean dose in dose areas per frame',obj.T.i));
        ylabel('mean dose')
        xlabel('frame number')
        grid on

        subplot(2,1,2)
        plot(dosemin)
        ylabel('minimum dose image value')
        xlabel('frame number')
        title(sprintf('if this has major dips, then some object \n rotates into the dose mask'))
        grid on
        fname=sprintf('Fig05 %02d dose area checks',obj.T.i);
        f.f_BoFig2PDF(h,fname)


    end

    %fit checks
    if obj.T.proofs==1
        h=figure(6);clf;
        imax=10;
        ii=ceil(obj.T.d.nframes(obj.T.cas,obj.T.rep)*rand(1,imax)); %check some random frames
        for i = 1:length(ii)

            im=squeeze(d(:,:,ii(i)));
            im(obj.T.mask==0)=[];
            rf =10; % reduction factor
            plot(obj.T.f{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
            hold on
            %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
            title(sprintf('%02d nun-uniformity data \n according to %d random frames',obj.T.i,imax ))
            xlabel('x')
            ylabel('y')
            zlabel('counts')
            pause(0.5)
        end

        fname=sprintf(sprintf('Fig06 %02d Beam fit',obj.T.i));
        f.f_BoFig2PDF(h,fname)

        %check if the corrected values are ok

        h=figure(7);clf;
        subplot(2,3,4:5);
        hold on
        ndf=size(obj.T.dose,2); % number of dose fields

        doseprofile=zeros(obj.T.d.nframes(obj.T.cas,obj.T.rep),ndf+1);
        oldprofile=zeros(obj.T.d.nframes(obj.T.cas,obj.T.rep),ndf+1);
        dosepix=zeros(1,ndf);

        for i=1:ndf %iterate dose fields
            % mean of each field for each frame, for both old (im) and
            % corrected (imc) data
            doseprofile(:,i)=squeeze(mean(mean((d_out(obj.T.dose(1,i):obj.T.dose(2,i),...
                obj.T.dose(3,i):obj.T.dose(4,i),:)))));
            oldprofile(:,i)=squeeze(mean(mean((d(obj.T.dose(1,i):obj.T.dose(2,i),...
                obj.T.dose(3,i):obj.T.dose(4,i),:)))));

            % plot
            plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))

            % the number of dose mask pixels
            dosepix(i)=(obj.T.dose(2,i)-obj.T.dose(1,i)+1)*(obj.T.dose(4,i)-obj.T.dose(3,i)+1);

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
        title(sprintf('%02d mean per frame dose of raw footage',obj.T.i))
        grid on



        subplot(2,3,[3,6]);
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        hold on
        for i=1:ndf
            r=rectangle('Position',[obj.T.dose(1,i),...
                obj.T.dose(3,i),...
                obj.T.dose(2,i)-obj.T.dose(1,i),...
                obj.T.dose(4,i)-obj.T.dose(3,i)],...
                'Facecolor','r');
        end
        title('dose areas')
        fname=sprintf('Fig07 %02d dose area checks',obj.T.i);
        f.f_BoFig2PDF(h,fname)
        obj.T.d.doseprofile{obj.T.i}=oldprofile(:,ndf+1);


    end
    T=rmfield(T,'unibeam');
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

    imshow(squeeze(im),[]);
    title(sprintf('Fig (%02d,%02d) Raw sinogram of the angle gauge region',cas,rep))
    fname=sprintf('Fig01 %02d %02d raw sinogram',cas,rep);
    %f.f_BoFig2PDF(h,fname)
end

function BoSave(d,basename,T)
    fname=sprintf('%s%s%02d.mat',obj.T.d.DataPath,basename,obj.T.i);
    fprintf('saving %s...',fname)
    savefast(fname,'d');
    fprintf('done. \n')
end

function BoSave2(array,fname,fpath)
    fprintf('saving %s...',fname)
    savefast(strcat(fpath,fname,'.mat'),'array');
    fprintf('done. \n')
end

function BoSave3(d,basename,T)
    fname=sprintf('%s%s%02d_%02d.mat',obj.T.d.DataPath,basename,obj.T.cas,obj.T.rep);
    fprintf('saving %s...',fname)
    savefast(fname,'d');
    fprintf('done. \n')
end

function BoSaveSino(d,T)
    %save function specifically for the data fom 2017 12 13
    jmax=[10,10,9,10];
    for i=1:4
        for j=1:jmax
            fname=fprintf('case_%02d_%02d.mat',i,j);
            BoSave2(d,obj.T.d.DataPath,fname)
        end
    end
end

function idiff=imdiff(im1,im2,x,y)
    a=abs(squeeze(im1-im2));
    idiff=sum(a(:));
end

function [shift,dif]=SinoMatch(ref,sino,T)
    disp(sprintf('matching case %02d',obj.T.i))
    % 1) do the rough finding via cross correlation
    guessF=f.SinoMatchXCOV(ref,sino,obj.T.match.GuessLine);

    % 2) do the fine matching iteratively
    [shift,dif]=f.SinoMatchIterate(ref,sino,0,guessF,T);
    
    

end

function T=FindShiftDiffArea(T)
    %% manually find the shift diff area
    % just plots the mask
    h=figure(106);clf

    
    DiffAreaMask=zeros(size(obj.T.sino.Raw(:,:,1)));
    imshow(squeeze(obj.T.sino.Raw(:,:,1)),[])
    DiffAreaMask([10:190,460:630],50:1450)=1;
    imshow(0.5*DiffAreaMask+obj.T.sino.Raw(:,:,1),[])
    obj.T.DiffAreaMask=DiffAreaMask;
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
    % f.test_SinoMatchXCOV(obj.T.sino.Raw(:,:,1),obj.T.sino.Raw(:,:,11),200)

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
    s1=obj.T.sino.Raw(:,:,1);
    col=hsv(39);
    for j=2:39 %iterate data
        s2=obj.T.sino.Raw(:,:,j);
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
    sxb=linspace(-5,5,21);  % basis iteration grd
    sfb=linspace(-5,5,21);  % i is i-dimenion, f is frames-dimension
    sx=sxb+guessX;                 % actual iteration grid
    sf=sfb+guessF;
    h=figure(10);
    iterations=3;
    for j=1:iterations % number of iterations
        diff2d=zeros(length(sx),length(sf)); % preallocate map
        for k=1:length(sx) % x dimension
            for l=1:length(sf) % frame dimension
                m=l+(k-1)*length(sx);
                %f.f_BoCount(m,20,10,5)
                im=f.fraccircshift(sino,[sx(k),sf(l)]);
                
                dif=abs(ref-im).*obj.T.DiffAreaMask;
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
            titlestring=sprintf('case %d \n %s',obj.T.i,titlestring);
        end
        title(titlestring)
        axis equal
        axis tight
        shift=[sx(I_row),sf(I_col)];
        dif=diff2d(I_row,I_col);
        
        % zoom in the shifts
        sxb=sxb/10;
        sfb=sfb/10;
        sx=sxb+sx(I_row);
        sf=sfb+sf(I_col);
    
    end
    fname=sprintf('Fig10 %02d sino shift match',obj.T.i);
    f.f_BoFig2PDF(h,fname)
    
    cminmax=[-0.1,0.1];
    h=figure(11);clf;
    subplot(2,1,1)
    im=f.fracshift(sino,shift);
    imshow(ref-im,[]);
    dif=abs(ref-im).*obj.T.DiffAreaMask;
    dif=sum(dif(:));
    title(sprintf('case %d difference after shifted (%.2f,%.2f), diff %.f',...
        obj.T.i,shift(1),shift(2),dif))
    colorbar
    caxis(cminmax)
    subplot(2,1,2)
    im2=f.fracshift(sino,[0,shift(2)]);
    imshow(ref-im2,[]);
    dif2=abs(ref-im2).*obj.T.DiffAreaMask;
    dif2=sum(dif2(:));
    title(sprintf('without x-shift, dif %.f',dif2))
    colorbar
    caxis(cminmax)
    fname=sprintf('Fig11 %02d sinogram match',obj.T.i);
    f.f_BoFig2PDF(h,fname)
    
    %M.shift(i,:)=[]
end

function [shift,mindif]=SinoMatchY(img,ref,guess,range,steps,mask,yscase)
    % special shift to correct for an y-offset in HS17 Xray Tomo
    % run test-yshift to check upon this routine
    % img= image to be shifted
    % ref= reference image
    % guess = manual guess for the shift in pixel
    % range = pixel range to be checked
    % steps = nr. of steps of the iterations
    % mask = difference evaluation area mask
    % yscase = case number for the y shift
    
    h=figure(12);clf;
    sb=linspace(-range,range,steps);
    s=sb+guess;
    iterations=3;
    for j=1:iterations % number of iterations
        dif2=zeros(length(s),1); % preallocate map
        for l=1:length(s) % frame dimension
            im=f.fraccircshift(img,[s(l),0]);
            dif=abs(ref-im).*mask;
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
            titlestring=sprintf('case %d \n %s',yscase,titlestring);
        end
        title(titlestring)
        % zoom in the shifts
        sb=sb/10;
        s=sb+shift;
    end
    
    
    fname=sprintf('Fig12 %02d y-shift correction',yscase);
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

function nframes=frames360(sino,refrow)
    % returns the number of frames in a sinogram that results
    % in the 360° cut
    s=size(sino);
    l=squeeze(sino(:,refrow));
    for i=1:s(2)
        diff(i)=sum((l-squeeze(sino(:,i))).^2);
    end
    diff(diff==0)=max(diff);
    [mindif,mini]=min(diff);
    nframes=mini-refrow;
    
    figure(111)
    plot(diff)
    title(sprintf('finding the 360° frames\n nframes=%d',nframes))
    xlabel('frame number')
    ylabel('difference')
    grid on
    
    h=figure(112);
    ax1=subplot(3,1,1);
    imshow(circshift(sino(:,refrow:refrow+nframes-2),[0,50]))
    title('-1 frame')
    ax2=subplot(3,1,2);
    imshow(circshift(sino(:,refrow:refrow+nframes-1),[0,50])) 
    title('correct cut')
    ax3=subplot(3,1,3);
    imshow(circshift(sino(:,refrow:refrow+nframes-0),[0,50]))    
    title('+1 frame')
    linkaxes([ax1,ax2,ax3])
    xlim([20,80])
    ylim([380,430])
    
    fname=('360 degree crop');
    f.f_BoFig2PDF(h,fname)
end

function [fits,fitshift,centershift]=FindCentering(T,d2,i)
    % finds the per plane shift required for centering
    % and fits a linear function over the channel height
    % T = T-struct
    % d2 = raw data block, such that d2(:,xx,:) is a sinogram
    % i = case iterator
    % fits = the parameters of the fit
    % fitsshift = the fitted shifts
    % centershift = the calculated shifts
    
    i=obj.T.i;
    k180=ceil(obj.T.q360.nFrames/2);
    centershift=zeros(1,1024);
    
    for j = 1:1024 % iterate planes
        % determine centering
        [c,lags]=xcov(d2(:,j,1),flipud(d2(:,j,k180)),40,'coeff');
        if any(isnan(c)) %debug
            %disp(sprintf('%4d contains NaN',j))
            continue
        else
            %disp(sprintf('%4d does not contain NaN',j))
            %ind_nan(k)=1;
        end
        [~,ind]=max(c);
        centershift(j)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
            (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
    end
    if i > 1
        fits=polyfit(400:1010,centershift(400:1010),1);
    else
        fits=polyfit(400:910,centershift(400:910),1); %case ref special
    end
    fitshift=polyval(fits,1:1024);
    %obj.T.c.centershift{i}=centershift;
    if obj.T.proofs==1
        f.CheckCenterFit(d2,fits,fitshift,centershift,i);
    end
    
end

function CheckCenterFit(d,fits,fitshift,centershift,i)
    h=figure(112);clf;
    plot(centershift,'Displayname', 'per-plane shift')
    hold on
    plot(fitshift,'Displayname','fitted shift')
    xlim([1,1024])
    ylim([-10,5])
    title(sprintf('%02d centering per y-plane and linear fit',i))
    xlabel('y plane (pixel planes)')
    ylabel('centering shift in pixel')
    grid on
    legend()
    fname=sprintf('Fig15 %02d centering shift',i);
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

function nframes=FindNframes(T,i)
    % returns the number of frames of .seq a video file
    fpath=strcat(obj.T.d.DataPath,obj.T.d.List(i).name);
    FileInfo=dir(fpath);
    FileSize=FileInfo.bytes;
    nframes=(FileSize-obj.T.d.header)/(obj.T.imsize(1)*obj.T.imsize(2)*obj.T.d.GreyBytes);
    if mod(nframes,1) ~= 0
        error('file size is not coherent with pixel sizes and header');
    end
end

function maxframes=FindMaxFrames(T)
    % this function does not work because the max frames number is
    % assessed only after sinograms of all data were made
    nframes=zeros(1,size(obj.T.d.List,1));
    for i = 1:size(obj.T.d.List,1)
        nframes(i)=f.FindNframes(T,i);
    end
    maxframes=max(nframes);
end

function y = loadSingleVariableMATFile(filename)
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
%  See also FIT, CFIT, SFIobj.T.

%  Auto-generated by MATLAB on 21-Jan-2018 12:20:44


%% Fit: 'untitled fit 1'.
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

function RemoveCameraSpots(d)
    %creates a mask that corrects the camera specific spots
    % manual spot findgin
    
    % 1) two dark spots x&y, box & center

    dscxy = [125,236;... % center x
             473,567;]  % center y
    mask=~mt.KernelGenBin(11,5,0);
     figure( 'Name', 'untitled fit 1' );
    for spot=1:size(dscxy,2)
        % create ringed data
        areax=(dscxy(1,spot)-5):(dscxy(1,spot)+5);
        areay=(dscxy(2,spot)-5):(dscxy(2,spot)+5);
        [x,y]=meshgrid(areax,areay);
        data=d(areax,areay);
        xr=x(mask);
        yr=y(mask);
        datar=data(mask)
         
        [fitresult{spot}, gof{spot}] = f.FitPixel1(xr, yr, datar)
        subplot(2,1,spot)
        h = plot( fitresult{spot}, [xr, yr], datar );
    
    end
    
   
   
    %legend( h, 'untitled fit 1', 'badimring vs. xring, yring', 'Location', 'NorthEast' );
    % Label axes
    xlabel xring
    ylabel yring
    zlabel data
    grid on
    view( -79.9, 11.6 );
    
end

function [Raw,bot,up,rad,yplane]=PreRead(T)
    toc1=toc;
    % reads in all data, and provides sinograms & radiograms
    % for manual parameter adjustment
    s=size(obj.T.d.List);
    is=obj.T.d.imsize;
    sy=obj.T.Raw.sinoyplanes;
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
            if obj.T.d.nframes(cas,rep) ~=0
                ind=rep+(cas-1)*maxrep;
                %disp([cas,rep,ind])
                d=f.ReadXrayTomo(T,cas,rep);
                
                % keep a copy of sinograms
                Raw{cas,rep}=squeeze(d(:,sy(1),:)); % corrected raw sinogram with angular marker
                bot{cas,rep}=squeeze(d(:,sy(2),:)); % bottom sinogram, above marker
                up{cas,rep}=squeeze(d(:,sy(3),:)); % high up sinogram
                rad{cas,rep}=squeeze(d(:,:,10)); %radiogram
                yplane{cas,rep}=squeeze(d(obj.T.Raw.yshiftplane,:,:));
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
    % finds start and end of the rotation
    % logic: first and last time the minimum difference between two
    % consecutive frames is larger than a threshold

    a=sino(:,2:end)-sino(:,1:end-1);
    mi=min(a,[],1);
    
    % find the flanks 
    start=find(abs(mi)>obj.T.Raw.CropThrash,1,'first')-obj.T.Raw.CropMargin;
    stop=find(abs(mi)>obj.T.Raw.CropThrash,1,'last')+obj.T.Raw.CropMargin;
    croprange=[start stop];
end

function CheckSino1(T,d)
if obj.T.proofs==1
        fid=1;
        h=figure(fid);clf;
        % raw sinogram plot
        subplot(2,1,1)
        imshow(squeeze(d(:,obj.T.Raw.sinoyplanes(1),:)),[]);
        ax=gca;
        hold on
        for i=1:2
            line(obj.T.Raw.CropRange(1,1,i)*[1,1],...
                [1,obj.T.d.imsize(2)],'color',[1,0,0])
        end
        tstr1=sprintf('(%02d,%02d) Raw sinogram of the angle gauge region',obj.T.cas,obj.T.rep);
        tstr2='with automated crop indicated; cropped & corrected sinogram';
        title(sprintf('%s \n %s',tstr1,tstr2));
%         xlabel('frames')
        %set(ax,'FontSize',14)
        %set(ax,'FontSize',14)
    %     tightInset = get(ax, 'TightInset');
%         position(1) = tightInset(1);
%         position(2) = tightInset(2);
%         position(3) = 1 - 2.6*tightInset(1) - tightInset(3);
%         position(4) = 1 - tightInset(2) - tightInset(4);
%         set(ax, 'Position', position);
%         set(ax,'units','centimeters')
%         pos = get(ax,'Position');
%         ti = get(ax,'TightInset');
%         
%         set(h, 'PaperUnits','centimeters');
%         set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%         set(h, 'PaperPositionMode', 'manual');
%         set(h, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%         
%         h.Renderer='Painters';

        
    end 
end

function CheckSino2(T,d)
     if obj.T.proofs==1
        fid=1;
        h=figure(fid);

        subplot(2,1,2);
        imshow(squeeze(d(:,obj.T.Raw.sinoyplanes(1),:)),[]);
        ax=gca;

        
        fname=sprintf('Fig01 %02d %02d raw sinogram',obj.T.cas,obj.T.rep);
        %set(ax,'FontSize',14)
    %     tightInset = get(ax, 'TightInset');
%         position(1) = tightInset(1);
%         position(2) = tightInset(2);
%         position(3) = 1 - 2.6*tightInset(1) - tightInset(3);
%         position(4) = 1 - tightInset(2) - tightInset(4);
%         set(ax, 'Position', position);
%         set(ax,'units','centimeters')
%         pos = get(ax,'Position');
%         ti = get(ax,'TightInset');
%         
%         set(h, 'PaperUnits','centimeters');
%         set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%         set(h, 'PaperPositionMode', 'manual');
%         set(h, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%         
%         h.Renderer='Painters';
        
        print(h,fname,'-dpdf')
        print(h,fname,'-dpng')
        savefig(h,fname)
        f.f_BoFig2PDF(h,fname)
        
    end 
end

function [d,sino,BS,MaskSino,MaskY,nframes]=SinoCrop11(T,d)
    % special because the cropped sino length is only know after the (1,1)
    % case. 
    
    % crop the first block!
    d=d(:,:,obj.T.Raw.CropRange(obj.T.cas,obj.T.rep,1):obj.T.Raw.CropRange(obj.T.cas,obj.T.rep,2));
    BS=size(d);
    %creates the Shiftmask xy
    MaskSino=f.MakeMask(obj.T.match.SinoWindows,...
        [obj.T.d.imsize(1),BS(3)]);
    MaskSino=f.MakeMask(obj.T.match.YWindows,...
        [obj.T.d.imsize(2),BS(3)]);
    %initialize sino array
    sino=zeros(size(d,1),size(d,3),obj.T.d.ncas,obj.T.d.nrep,size(obj.T.Raw.sinoyplanes,2),'single');
    
    %copy sinograms
    for i = 1:2
        sino(:,:,obj.T.cas,obj.T.rep,i)=squeeze(d(:,obj.T.Raw.sinoyplanes(i),:));
    end
    nframes=BS(3);
end

function [d,sino,nframes]=SinoCrop(T)
    % adjust to (1,1) sinogram length
    
    d=d(:,:,obj.T.Raw.CropRange(cas,rep,1):obj.T.Raw.CropRange(cas,rep,2));
    sino=zeros(size(d,1),size(d,3),1,1,size(obj.T.Raw.sinoyplanes));
    for i = obj.T.Raw.sinoyplanes
        sino(:,:,1,1,i)=squeeze(d(:,obj.T.Raw.sinoyplanes(i),:));
    end
    nframes=size(d,3);
end

function d=SinoLengthAdjust(d,T)
    % adjust the 3rd dimension of a sino to that of a reference sino. 
    % d the data block
    % T the T struct
    
    % 1) get length of the d after realcrop
    rc=squeeze(obj.T.Raw.CropRange(obj.T.cas,obj.T.rep,:));
    dLen=rc(2)-rc(1)+1;
    diff=size(obj.T.Raw.Sino,2)-dLen;
    
    if mod(diff,1)~=0 % its not an integer!
        error('sino length diff is not integer! rc=(%.2f,%.2f), diff=%.2f',rc(1),rc(2),diff)
    end

    a=ceil(diff/2);
    b=floor(diff/2);
    rc2=rc+[-a,b]; 
    
    if rc2(1)<1 % need to pad ath the beginning
        % pad by the difference
        padLen=1-rc2(1);
        d=cat(3,repmat(d(:,:,1),[1,1,padLen]),d);
        %correct rc
        rc2(1)=1;
    end
        
    if rc2(2)>size(d,3) % need to pad at the end
        % pad
        padLen=rc2(2)-size(d,3);
        d=cat(3,d,repmat(d(:,:,end),[1,1,padLen]));
        % correct rc
        rc2(2)=size(d,3);
    end
    
    % finally crop!
    d=d(:,:,rc2(1):rc2(2));

    
end

function [d_out,fi]=BeamAndDoseCorrection2(d,T)
    % corrects for beam & dose
    % d = data block
    % T = parameters struct
    fprintf('%d %d correctiong for beam and dose...\n',obj.T.cas,obj.T.rep)
    mask=obj.T.dose.mask; %re-label
    
    % pick the valid points from the mask & prepare for fitting
    maskmesh=size(mask);
    [x,y]=ndgrid(1:maskmesh(1),1:maskmesh(2)); %make grid
    xx=x;
    yy=y;
    xx(mask==0)=[]; % remove pixels outside of dose areas
    yy(mask==0)=[];

    unibeam=ones(size(d),'single'); % contains the correction factor

    % find the fits
    %fprintf('%d %d finding beam non-uniformity fits...\n',obj.T.cas,obj.T.rep)
    nframes=obj.T.Raw.nframes(obj.T.cas,obj.T.rep); % parfor-hack
    parfor j=1:nframes

        f.f_BoCount(j,50,10,5)
        im=squeeze(d(:,:,j));    % check out a slice from stack
        im(mask==0)=[];             % remove non-flatfield-areas
        fi{j}= fit( [xx',yy'],double(im'), 'poly22' );   % make a fit
        unibeam(:,:,j)=feval(fi{j},x,y);         % make the beam

    end
   
    fprintf('\n')

    % correct the non uniformity
    d_out=d./unibeam;

    % Dose Area checks
    if obj.T.proofs==1 % check the areas

        h=figure(4);clf;
        set(gcf,'name','dose area mask check')
        subplot(1,3,1)
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d original',obj.T.cas,obj.T.rep))
        subplot(1,3,2)
        imshow(obj.T.dose.mask',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d mask',obj.T.cas,obj.T.rep))
        subplot(1,3,3)
        imshow((single(squeeze(d(:,:,1))).*(0.5+obj.T.dose.mask/2))',[]);set(gca,'YDir','normal')
        title(sprintf('%02d %02d combined',obj.T.cas,obj.T.rep))
        fname=sprintf('Fig04 %02d dose mask',obj.T.cas,obj.T.rep);
        f.f_BoFig2PDF(h,fname)

        %check if nothing rotates in these areas
        h=figure(5);clf;
        dosemin=zeros(1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep));
        dosemean=zeros(1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep));
        for i=1:obj.T.Raw.nframes(obj.T.cas,obj.T.rep)
            a=single(squeeze(d(:,:,i))).*obj.T.dose.mask;
            dosemin(i)=min(a(a>0));
            dosemean(i)=mean(a(a>0));
        end

        subplot(2,1,1)
        plot(dosemean)
        title(sprintf('%02d %02 mean dose in dose areas per frame',obj.T.cas,obj.T.rep));
        ylabel('mean dose')
        xlabel('frame number')
        xlim([1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep)])
        grid on

        subplot(2,1,2)
        plot(dosemin)
        ylabel('minimum dose image value')
        xlabel('frame number')
        title(sprintf('if this has major dips, then some object \n rotates into the dose mask'))
        grid on
        xlim([1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep)])
        fname=sprintf('Fig05 %02d %02d dose area checks',obj.T.cas,obj.T.rep);
        f.f_BoFig2PDF(h,fname)


    end

    %fit checks
    if obj.T.proofs==1
        h=figure(6);clf;
        imax=1;
        ii=ceil(obj.T.Raw.nframes(obj.T.cas,obj.T.rep)*rand(1,imax)); %check some random frames
        for i = 1:length(ii)

            im=squeeze(d(:,:,ii(i)));
            im(obj.T.dose.mask==0)=[];
            rf =10; % reduction factor
            plot(fi{ii(i)},[xx(1:rf:end)',yy(1:rf:end)'],im(1:rf:end)')
            hold on
            %plot(f{i},[x(1:100:end),y(1:100:end)],im(1:100:end))
            title(sprintf('%02d %02d nun-uniformity data \n according to a random frame',obj.T.cas,obj.T.rep))
            xlabel('x')
            ylabel('y')
            zlabel('counts')
            %pause(0.5)
        end

        fname=sprintf(sprintf('Fig06 %02d %02d Beam fit',obj.T.cas,obj.T.rep));
        f.f_BoFig2PDF(h,fname)

        %check if the corrected values are ok

        h=figure(7);clf;
        subplot(2,3,4:5);
        hold on
        ndf=size(obj.T.dose.windows,2); % number of dose fields

        doseprofile=zeros(obj.T.Raw.nframes(obj.T.cas,obj.T.rep),ndf+1);
        oldprofile=zeros(obj.T.Raw.nframes(obj.T.cas,obj.T.rep),ndf+1);
        dosepix=zeros(1,ndf);

        for i=1:ndf %iterate dose fields
            % mean of each field for each frame, for both old (im) and
            % corrected (imc) data
            doseprofile(:,i)=squeeze(mean(mean((d_out(...
                obj.T.dose.windows(1,i):obj.T.dose.windows(2,i),...
                obj.T.dose.windows(3,i):obj.T.dose.windows(4,i),:)))));
            oldprofile(:,i)=squeeze(mean(mean((d(...
                obj.T.dose.windows(1,i):obj.T.dose.windows(2,i),...
                obj.T.dose.windows(3,i):obj.T.dose.windows(4,i),:)))));

            % plot
            plot(doseprofile(:,i),'DisplayName',sprintf('region %d',i))

            % the number of dose mask pixels
            dosepix(i)=(obj.T.dose.windows(2,i)-obj.T.dose.windows(1,i)+1)...
                *(obj.T.dose.windows(4,i)-obj.T.dose.windows(3,i)+1);

            % the total dose, should be constantly equal 1
            doseprofile(:,ndf+1)=doseprofile(:,ndf+1)+dosepix(i)*doseprofile(:,i);
            oldprofile(:,ndf+1)=oldprofile(:,ndf+1)+dosepix(i)*oldprofile(:,i);
        end
        xlim([1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep)])
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
        title(sprintf('%02d %02d mean per frame dose of raw footage',obj.T.cas,obj.T.rep))
        grid on
        xlim([1,obj.T.Raw.nframes(obj.T.cas,obj.T.rep)])

        subplot(2,3,[3,6]);
        imshow(squeeze(d(:,:,1))',[]);set(gca,'YDir','normal')
        hold on
        for i=1:ndf
            r=rectangle('Position',[obj.T.dose.windows(1,i),...
                obj.T.dose.windows(3,i),...
                obj.T.dose.windows(2,i)-obj.T.dose.windows(1,i),...
                obj.T.dose.windows(4,i)-obj.T.dose.windows(3,i)],...
                'Facecolor','r');
        end
        title('dose areas')
        fname=sprintf('Fig07 %02d %02d dose area checks',obj.T.cas,obj.T.rep);
        f.f_BoFig2PDF(h,fname)
        


    end
   
end

function [rad,croprange,nframes,shiftXZ,diffXZ]=PreProcessingPreAllocation(T)
    % returns empty arrays. 
    % main pupose is to keep the main script tidy
    rad=zeros([obj.T.d.imsize,obj.T.d.ncas,obj.T.d.nrep],'single');    % radiography
    croprange=zeros(obj.T.d.ncas,obj.T.d.nrep,2); % stores the crop range data
    nframes=zeros(obj.T.d.ncas,obj.T.d.nrep); %nframes after cropping
    shiftXZ=zeros(obj.T.d.ncas,obj.T.d.nrep,3); % XYZ shifts from matching
    diffXZ(obj.T.d.ncas,obj.T.d.nrep);  % XZ shift differences
    
end

function CheckSinoShiftMask(T)
    %% manually find the shift diff area
    % just plots the mask
    h=figure(106);clf
    im=squeeze(obj.T.Raw.Sino(:,:,1,1,1)).*(0.5*obj.T.match.DiffMaskSino+0.5);
    imshow(im,[])
    fname=sprintf('shift base area plot');
    f.f_BoFig2PDF(h,fname)
end

function s=TimeString()
    % returns a string of the current time
    s=string()
    
end

end %static
end %class










% publication ffigure stuff

% %%
%     set(ax,'FontSize',14)
%     tightInset = get(ax, 'TightInset');
%     position(1) = tightInset(1);
%     position(2) = tightInset(2);
%     position(3) = 1 - 2.6*tightInset(1) - tightInset(3);
%     position(4) = 1 - tightInset(2) - tightInset(4);
%     set(ax, 'Position', position);
%     set(ax,'units','centimeters')
%     pos = get(ax,'Position');
%     ti = get(ax,'TightInset');
%     
%     set(h, 'PaperUnits','centimeters');
%     set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%     set(h, 'PaperPositionMode', 'manual');
%     set(h, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%     
%     h.Renderer='Painters';
%     
%     print(h,fname,'-dpdf')
%     print(h,fname,'-dpng')
%     savefig(h,fname)




















