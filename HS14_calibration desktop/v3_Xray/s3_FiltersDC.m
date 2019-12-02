
%% DC filtering and adding
%check if we shall run this script
if  any(itr_stp==3) % execute this script?
    i=3; % specifies step index of this file, is two
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile
    
    %% manual spot filter for DC
    
    for j = 4
        % iterate emp, d2o, cl3, dc
        log_wrk=logfile{i-1}{j}; %working copy of the logfile
        n_blocks=max(log_wrk{:,7}); % number of blocks for this case
        
        disp(strcat(nam_block,nam_stp{i},'_',pfx_cas{j}));
        %load the block from previous step
        block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{j});
        path_block=strcat(dir_work,block_old);
        block=f3_loadSingleVariableMATFile(path_block);
        [x,y,z]=size(block);


        for k=1:z
            im=block(:,:,k);
            %filter for bright & dark spots


            block(:,:,k)=im;
        end %k-loop ends

        DC=sum(block,3)/z;

        %save block to disk
        if ind_write==1
            savefast(strcat(dir_work,nam_DC),'DC');
        end

            

        
    end % j-loop ends

    % j ends
    clear block    
       % load
       % open file
       % append block
    
end
