%% Preprocessing 2: only modify the shifts
% delete(gcp('nocreate'));
% parpool(7);
% preprocess
tic
for cas=[7,8]%]1:T.d.ncas %iterate cases
    if any(T.d.nframes(cas,:)~=0) % skip if no files available
    T.cas=cas; % copy iterater to T struct
    d_add=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
    d_add2=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
    for rep=1:T.d.nrep % iterate repetitions
        d=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
        d2=zeros([T.Raw.BS(1:2),T.q360.nFrames(1,1)],'single');
        if T.d.nframes(cas,rep)~=0 % skip if no measurement available
            T.rep=rep; % copy iterater to T struct
            fprintf('%s shift- modifying\n',f.CommLineStart(cas,rep)); %debug

            
            % load intermediate result
            d=f.BoLoad(T.fnames.corr{cas,rep},T);
            
            % shift corrections
            if all([T.cas,T.rep]==[1,1]) % the first [1,1]
                % first sinogramm, nothing to shift here
                T.Match.YPlanes=zeros(T.Raw.BS(2),T.Raw.BS(3),...
                    T.d.ncas,T.d.nrep,T.sys.StoreFormat);
                d2=d;
            end
            T.Match.YPlanes(:,:,cas,rep)=...
                squeeze(d(T.Match.yshiftplane,:,:));
            
            if ~all([T.cas,T.rep]==[1,1]) % all other cases
                % shift towards (1,1)in x & z direction
                [T.Match.shift(cas,rep,1:2),T.Match.diff(cas,rep)]=...
                    f.SinoMatch(squeeze(T.Raw.Sino(:,:,1,1,1)),...
                    squeeze(d(:,T.Raw.sinoyplanes(1),:)),T);
                
%                %apply x and z shift
% old             d=f.fraccircshift(d,[T.Match.shift(cas,rep,1),...
% old                 0,T.Match.shift(cas,rep,2)]);
               % apply only z shift
                 d=f.fraccircshift(d,[T.Match.shift(cas,rep,1),...
                    0,T.Match.shift(cas,rep,2)]);
                d2=f.fraccircshift(d,[0,...
                    0,T.Match.shift(cas,rep,2)]);
                
%                 % yshift
%                 [T.Match.shift(cas,rep,3),T.Match.diff(cas,rep,3)]=...
%                     f.SinoMatchY(T);
%                 
%                 % apply shift
%                 d=f.fraccircshift(d,[0,T.Match.shift(cas,rep,3),0]);
                
            end
            
            % crop z
            d=d(:,:,T.q360.startframe:(...
                T.q360.startframe+T.q360.nFrames(1,1)-1));
            d2=d2(:,:,T.q360.startframe:(...
                T.q360.startframe+T.q360.nFrames(1,1)-1));
            
            
        end
        d_add=d_add+d; % sum up the data blocks
        d_add2=d_add2+d2; % sum up the data blocks
        close all
        
        
    end %reps
    clear d d2 % free up memory
    d_add=d_add/(sum(T.d.nframes(cas,:)~=0)); % normalize back ot [0 1]
    d_add2=d_add2/(sum(T.d.nframes(cas,:)~=0)); % normalize back ot [0 1]
   
    T.fnames.addWS{cas}=f.BoSave4(d_add,'2_addWS_',T); % save file & filename
    T.fnames.addWOS{cas}=f.BoSave4(d_add2,'2_addWOS_',T); % save file & filename
    
    if cas~=1 %

        d_add1 = f.BoLoad(T.fnames.add{1},T); %load ref
        d_add=-log(d_add./d_add1); % Beer-Lambert referencing
        d_add2=-log(d_add2./d_add1); % Beer-Lambert referencing
        T.Log.Sino(:,:,cas)=d_add(:,T.Raw.sinoyplanes(1),:); % keep a real sino
        clear d_add1
        T.fnames.divWS{cas}=f.BoSave4(d_add,'3_divWS_',T);
        T.fnames.divWOS{cas}=f.BoSave4(d_add2,'3_divWOS_',T);
    end
    clear d_add d_add2
    fprintf('%s saving T...',f.CommLineStart2(cas)); %debug
    save('T','T') %save T struct
    fprintf('done. \n')
    
    end
end %cases % preprocessing 2

%% centering & reconstruction
% find centering

T.Cen.start=173; % rough cut to center, in terms of x-pixel
T.Cen.stop=491;
T.Cen.range=T.Cen.start:T.Cen.stop; % in here lies the channel
T.Rec.recsize=length(T.Cen.range); % reconstruction edge size
T.Rec.angles=linspace(0,360,T.q360.nFrames(T.cas,T.rep)+1); %make the angles list
T.Rec.angles(end)=[];
T.Cen.CSFR{2}=345:1010; % CentShiftFitRange: y-planes, based on which the 
T.Cen.CSFR{1}=400:910;  % center shift will be performed
T.Rec.SinoRotShift=200; % sino is circshifted by this to have the 4 quadrants
                        % aligned

% delete(gcp('nocreate'));
%% parpool(4);
k=1;
recontot=[];
for cas =[1,6,7,8]%:T.d.ncas
if ~isempty(T.fnames.add{cas}) % skip if no files available
    T.cas=cas;
    
    % preallocate the recon block
    recon=zeros(T.Rec.recsize,T.Rec.recsize,T.d.imsize(2),'single'); % block of recons
    if cas>1 % usual case, take div files
        rfnameW=T.fnames.divWS{cas}; 
        rfnameWO=T.fnames.divWOS{cas}; 
        %rfname=T.fnames.divNewShift{cas}; 
        
    else % for case 1, reconstruct the absolute file
        rfnameW=T.fnames.addWS{cas}; 
         rfnameWO=T.fnames.addWOS{cas}; 
        %rfname=T.fnames.addNewshift{cas};
    end
    
    % load sinograms
    d=f.BoLoad(rfnameW,T);
    d2=f.BoLoad(rfnameWO,T);
    
    % crop to specified range
    d=d(T.Cen.range,:,:);
    d2=d2(T.Cen.range,:,:);

    % find the centering for each frame
    [fit1{cas},fitshift1(cas,:),centshift1(cas,:)]=...
        f.FindCentering(T,d); 
    f2=copyobj(gcf,0);
    gca()
    title(sprintf('%s\n%s',f.FigTName2(...
        'centering WS',13,T.cas),...
        sprintf('shift from %.1f to %.1f, diff %.3f',...
        fitshift1(1),fitshift1(end),fitshift1(end)-fitshift1(1))))
    fname=f.FigFileName2('centering shiftWS',13,T.cas);
    f.f_BoFig2PDF(f2,fname)
    
    
    [fit2{cas},fitshift2(cas,:),centshift2(cas,:)]=...
        f.FindCentering(T,d2); 
    title(sprintf('%s\n%s',f.FigTName2(...
        'centering WOS',13,T.cas),...
        sprintf('shift from %.1f to %.1f, diff %.3f',...
        fitshift2(1),fitshift2(end),fitshift2(end)-fitshift2(1))))
    f3=copyobj(gcf,0);
    fname=f.FigFileName2('centering shiftWOS',13,T.cas);
    f.f_BoFig2PDF(f3,fname)
    
    % turn off warning from iradon about NaNs in lines
    id='MATLAB:interp1:NaNstrip'; warning('off',id);
    
    % parforhacks
    sinsh=T.Rec.SinoRotShift;
    fsh1=fitshift1;
    fsh2=fitshift2;
    ang=T.Rec.angles;
    recsize=T.Rec.recsize;
    BS=T.Raw.BS;
    fprintf('planes: ')
    
    for plane=[50,500,1000]%1:BS(2)
        %disp(plane)
        f.f_BoCount(plane,20,10,5)
        try
            %rotate the data such that the rods are in the quadrants
            imc1=circshift(squeeze(d(:,plane,:)),sinsh,2); 
            imc2=circshift(squeeze(d2(:,plane,:)),sinsh,2); 
            
            
            % center
            imc=f.fraccircshift(imc,-fsh1(cas,plane));
            imc2=f.fraccircshift(imc2,-fsh2(cas,plane));
            
            % reconstruct
            rec1=iradon(imc1,ang,'spline','Hann',1,recsize);
            rec2=iradon(imc2,ang,'spline','Hann',1,recsize);
            %disp(sprintf('%d %d',j,size(imc,1)));
            recon(:,:,plane)=rec1; %hardcode after-recon-cropping
            recontot1(k,:,:)=rec1;
            recontot2(k,:,:)=rec2;
            k=k+1;
        catch e %e is an MException struct
            fprintf(2,'%s\n',e.identifier);
            fprintf(2,'%s\n',e.message);
            fprintf('%4d didnt reconstruct\n',plane)
        end
    end
    T.fnames.rec{cas}=f.BoSave4(recon,'4_rec_',T);
       
end   % isempty(T.fnames.add{cas}) % skip if no files available 
end
%%


