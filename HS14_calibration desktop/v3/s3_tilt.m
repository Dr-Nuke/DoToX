%% s3_tilt for correcting th tilt

if  any(itr_stp==6) % execute this script?
    i=6; % specifies step index of this file, is two
    t1=toc;
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    disp('tilt correction...')
   
    % load suited data set for tilt detection 
    for j=[2,3] % only D2o and cl3

        block_old=strcat(nam_block,nam_stp{4},'_',pfx_cas{j});
        path_block=strcat(dir_work,block_old);
        tiltblock=f3_loadSingleVariableMATFile(path_block);

        for k=1:size(angles,2) % 375 für the HS14 tomo
            tilt(j,k)=f3_subpixcorr(tiltblock(:,tiltline(1),k),...
                                             tiltblock(:,tiltline(2),k));
        end
    end
    
    tiltheight=tiltline(2)-tiltline(1)+1;
    meantilt=atand(mean(tilt(:))/(tiltheight));

    %load and modify the 5div data
    block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{2});
    path_block=strcat(dir_work,block_old);
    tiltblock=f3_loadSingleVariableMATFile(path_block);    

    % add here the rotation
    % deactivated this for now. is not working well

    %tiltblock=imrotate(tiltblock,-meantilt,'bilinear','crop');


    % generate variable name
    v=strcat(nam_block,nam_stp{i},'_',pfx_cas{3});
    eval([v '= tiltblock;']);
    %save block to disk
    if ind_write==1
        fprintf('saving %s to harddrive.....',v)
        savefast(strcat(dir_work,v),v);
        fprintf('done\n')
    end
    eval(['clear ' v]); % delete from workspace



    %clear tiltblock
    t2=toc;
    fprintf(t_string,name,t2-t1,t2/60);

        
    
end

    
   