
T.Cen.start=175; % rough cut to center, in terms of x-pixel
T.Cen.stop=481;
T.Cen.range=T.Cen.start:T.Cen.stop; % in here lies the channel
T.Rec.recsize=length(T.Cen.range); % reconstruction edge size
T.Rec.angles=linspace(0,360,T.q360.nFrames(T.cas,T.rep)+1); %make the angles list
T.Rec.angles(end)=[];
T.Cen.CSFR{2}=345:1010; % CentShiftFitRange: y-planes, based on which the 
T.Cen.CSFR{1}=400:910;  % center shift will be performed
T.Rec.SinoRotShift=200; % sino is circshifted by this to have the 4 quadrants
                        % aligned


for cas =1
if ~isempty(T.fnames.add{cas}) % skip if no files available
    T.cas=cas;
    
    % preallocate the recon block
    recon=zeros(T.Rec.recsize,T.Rec.recsize,T.d.imsize(2),'single'); % block of recons
    if cas>1 % usual case, take div files
        rfname=T.fnames.div{cas}; 
    else % for case 1, reconstruct the absolute file
        rfname=T.fnames.add{cas}; 
    end
    
    % load sinograms
    d=f.BoLoad(rfname,T);
    
    % crop to specified range
    d=d(T.Cen.range,:,:);
    
    % find the centering for each frame
    [T.Cen.fit{cas},T.Cen.fitshift(cas,:),T.Cen.centshift(cas,:)]=...
        f.FindCentering(T,d); 
end
end