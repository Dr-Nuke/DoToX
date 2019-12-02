%% files adder
% adds those slices that were recorded multiple times


%%filtering
%check if we shall run this script
if  any(itr_stp==4) % execute this script?
    t1=tic
    i=4; % specifies step index of this file, is two
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile
    
    
    for j = [2,3]  % iterate  d2o, cl3, NOT dc and emp
        log_wrk=logfile{i-1}{j}; %working copy of the logfile
        n_blocks=max(log_wrk{:,7}); % number of blocks for this case
        block_added=zeros(imrange(1,2),imrange(2,2),nangles,'single'); %+1 for 180� image

        disp(strcat(nam_block,nam_stp{i},'_',pfx_cas{j}));
        %load the block from previous step
        block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{j});
        path_block=strcat(dir_work,block_old);
        block=f3_loadSingleVariableMATFile(path_block);
        [x,y,z]=size(block);
        
        
        zs=1:num_proj:z; % start slice for current projection
        zf=zs+(num_proj-1); % finish slice for current projection
        
        
        for k=1:length(zs)
            f3_BoCountDisp(k,300)
            projNr=log_wrk{zs(k),2};
            if all(log_wrk{zs(k):zf(k),2}==projNr) %check if all 
                % projection numbers are equal
                

                    %block_added(:,:,nangles)=sum(block(:,:,zs(nangles):zf(nangles)),3);

                    block_added(:,:,k)=sum(block(:,:,zs(k):zf(k)),3);
               
            end
            
        end %k-loop ends

        % generate variable name
        v=strcat(nam_block,nam_stp{i},'_',pfx_cas{j});
        eval([v '= block_added;']);

        %save block to disk
        if ind_write==1
            fprintf('saving %s to harddrive.....',v)
            savefast(strcat(dir_work,v),v);
            fprintf('done\n')
        end
        eval(['clear ' v]); % delete from workspace

    end % j-loop ends

    % j ends
    clear block block_added
       % load
       % open file
       % append block
    

    t2=toc;
    fprintf(t_string,name,t2-t1,t2/60);
end