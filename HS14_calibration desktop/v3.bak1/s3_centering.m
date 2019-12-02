%% s3_centering 

if  any(itr_stp==7) % execute this script?
    i=7; % specifies step index of this file, is two
    t1=toc;
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
        
    
       
    % cl3
    block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{3});
    path_block=strcat(dir_work,block_old);
    block=f3_loadSingleVariableMATFile(path_block);

    recsize=999;
    ind_nan=zeros(1,imrange(2,2));
    recontot=zeros(751,751,imrange(2,2),'single');
    
    centrange=[360:2090]; % the relevant range for centering in y pixel 
    shift=zeros(1,imrange(2,2));
 %%   


% determine centering 
    for k=1:2450 %1:imrange(2,2)
        % impix{2} %1:imrange(2,2) %for each pixel plane
        %disp(k)

        sino=squeeze(block(:,k,:)); %sinogram

        % find the cross correlation maximum
        [c,lags]=xcov(sino(:,1),flipud(sino(:,end)),'coeff');

        %check if centring works
        if ind_nan(k)~=0
            disp(sprintf('%4d already done before',k))
            %continue
        end
        
        if any(isnan(c))
            disp(sprintf('%4d contains NaN',k))
            ind_nan(k)=-1;
            continue
            
        else
            disp(sprintf('%4d does not contain NaN',k))
            ind_nan(k)=1;
        end


        [~,ind]=max(c); 



        % the three point gaussian maximum fit
        shift(k)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/(log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));

    end
    centfit=polyfit(centrange,shift(centrange),1);
    shift_orig=shift;
    shift=polyval(centfit,1:imrange(2,2));
    
%%

    % center correction
    for l=1:24
        a=(l-1)*100+1;
        b=l*100;
        disp([l,a,b])
        parfor k=a:b
        disp(k)    
            
        sino=squeeze(block(:,k,:)); %sinogram
        sino=f3_center(sino,shift(k));

        %hardcode cropping
        sino=sino(201:end-200,:);
        try
            recon=iradon(circshift(sino(:,1:end-1),11,2),angles,'spline','Hann',1,999);
            recontot(:,:,k)=recon(100:850,100:850); %hardcode after-recon-cropping
        catch
            disp(sprintf('%4d didnt reconstruct',k))
        end
        %sinotot(:,:,k)=sino;
        end





    end
    %%
    
     % generate variable name
    v=strcat(nam_block,nam_stp{i},'_',pfx_cas{3});
    eval([v '= recontot;']);
    %save block to disk
    if ind_write==1
        fprintf('saving %s to harddrive.....',v)
        savefast(strcat(dir_work,v),v);
        fprintf('done\n')
    end
    eval(['clear ' v]); % delete from workspace


    %clear block inte frac bs XI YI AI c lags   
    t2=toc;
    fprintf(t_string,name,t2-t1,t2/60);  
end

    
