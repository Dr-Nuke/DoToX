
%%filtering
%check if we shall run this script
if  any(itr_stp==3) % execute this script?
    t1=toc
    i=3; % specifies step index of this file, is two
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile
    
    %% manual spot filter 
    
    
    for j = [1,2,3]  % iterate emp, d2o, cl3, dc
        log_wrk=logfile{i-1}{j}; %working copy of the logfile
        n_blocks=max(log_wrk{:,7}); % number of blocks for this case

        disp(strcat(nam_block,nam_stp{i},'_',pfx_cas{j}));
        
        %load the block from previous step
        block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{j});
        path_block=strcat(dir_work,block_old);
        block=f3_loadSingleVariableMATFile(path_block);
        [x,y,z]=size(block);


        parfor k=1:z
            f3_BoCountDisp(k,300)
            im=block(:,:,k);
            %filter for bright & dark spots
            im=f_MediFilter2( im,kernel_medi2,thresh_medi2,'median' );

            % filter gammaspots(median)
            % after-filtering to remove circles left overy by spotfilter
            im=f_MediFilter1(im,thresh_medi1);    %vertical
            im=f_MediFilter1_2(im,thresh_medi1);  % horizontal (or vice versa) 

            block(:,:,k)=im;
        end %k-loop ends

        %generate variable name
        v=strcat(nam_block,nam_stp{i},'_',pfx_cas{j});
        eval([v '= block;']);

        
        if ind_write==1
            fprintf('saving %s to harddrive.....',v)
            savefast(strcat(dir_work,v),v);
            fprintf('done\n')
        end
        eval(['clear ' v]); % delete from workspace
            

    % write down log file
    logfile{i}{j}=log_wrk;
      
    end % j-loop ends
    %save logfile to disk
    save(nam_log_ful,'logfile');
 

    % clear block    
       % load
       % open file
       % append block
    t2=toc;
    fprintf(t_string,name,t2-t1,t2/60);
end
