%% block forming script


%check if we shall run this script
if  any(itr_stp==2) % execute this script?    
    t1=toc;
    i=2; % specifies step index of this file, is two
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile

    for j = itr_cas   % iterate emp, d2o, cl3, dc
        log_wrk=logfile{i-1}{j};
        log_wrk.Properties.VariableNames(7)={'blockNr'}; %reuse column
        log_wrk.Properties.VariableNames(8)={'BlockSliceNr'}; %reuse column
        
        k_tot=size(log_wrk,1); % number of files
        disp(strcat(pfx_cas{j},': ', num2str(k_tot),' files'));
        
        dir_r=strcat(dir_work,dir_stp{i-1},dir_cas{j}); % read folder
        dir_w=dir_work ;  % write folder
        
       
        % k is the image number
        
        block=zeros(imrange(1,2),imrange(2,2),k_tot,'single');
        
        parfor k=1:k_tot
            f3_BoCountDisp(k,300)
            file_r=log_wrk{k,14}{1};
            filename=strcat(dir_r,file_r);
            block(:,:,k)=im2single(fitsread(filename)');
        end
        
        %generate variable name
        v=strcat(nam_block,nam_stp{i},'_',pfx_cas{j});
        disp(strcat('writing ',v));
        eval([v '= block;']);

        %save block to disk
        if ind_write==1
            savefast(strcat(dir_work,v),v);
        end
    eval(['clear ' v]); % delete from workspace

    % write down log file
    logfile{i}{j}=log_wrk;
      
    end % j-loop ends
    %save logfile to disk
    save(nam_log_ful,'logfile');
    
    
    % j ends
    clear block    
    t2=toc;
    fprintf(t_string,name,t2-t1,t2/60)
end