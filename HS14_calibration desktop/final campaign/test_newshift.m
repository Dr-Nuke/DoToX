recona=zeros(T.Rec.recsize,T.Rec.recsize,T.d.imsize(2),T.d.ncas,'single');
%%
for cas =1%2:T.d.ncas
if ~isempty(T.fnames.add{cas}) % skip if no files available
    T.cas=cas;
    
    % preallocate the recon block
    reconcas=zeros(T.Rec.recsize,T.Rec.recsize,T.d.imsize(2),'single'); % block of recons

        rfname=T.fnames.add{cas}; 
        %rfname=T.fnames.addNewshift{cas};
    
    % load sinograms
    d=f.BoLoad(rfname,T);
    
    % crop to specified range
    d=d(T.Cen.range,:,:);
    
    % find the centering for each frame
    [T.Cen.fit{cas},T.Cen.fitshift(cas,:),T.Cen.centshift(cas,:)]=...
        f.FindCentering(T,d); 
    
    % turn off warning from iradon about NaNs in lines
    id='MATLAB:interp1:NaNstrip'; warning('off',id);
    
    % parforhacks
    sinsh=T.Rec.SinoRotShift;
    fsh=T.Cen.fitshift;
    ang=T.Rec.angles;
    recsize=T.Rec.recsize;
    BS=T.Raw.BS;
    fprintf('planes: ')
    
    for plane=[50,100,500,1000]%
        %disp(plane)
        f.f_BoCount(plane,20,10,5)
        try
            %rotate the data such that the rods are in the quadrants
            imc=circshift(squeeze(d(:,plane,:)),sinsh,2); 
            
            
            % center
%            imc=f.fraccircshift(imc,-fsh(cas,plane));
            imc=f.fraccircshift(imc,-fsh(7,plane));
            imc=-log(imc);
            
            % reconstruct
            rec=iradon(imc,ang,'spline','Hann',1,recsize);
            %disp(sprintf('%d %d',j,size(imc,1)));
            reconcas(:,:,plane)=rec; %hardcode after-recon-cropping
            %recontot(k,:,:)=rec;
            %k=k+1;
        catch e %e is an MException struct
            fprintf(2,'%s\n',e.identifier);
            fprintf(2,'%s\n',e.message);
            fprintf('%4d didnt reconstruct\n',plane)
        end
    end
    recona(:,:,:,cas)=reconcas;
    
       
end   % isempty(T.fnames.add{cas}) % skip if no files available 
end
%%

