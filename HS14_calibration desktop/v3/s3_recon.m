%% s3_recon


if  any(itr_stp==8) % execute this script?
    i=8; % specifies step index of this file, is eight
    t1=toc;
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile
    
    
    
    j=3; % cl3
    block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{j});
    path_block=strcat(dir_work,block_old);
    block=f3_loadSingleVariableMATFile(path_block);
    
    for k=1000 %impix{2}
        
        recon=iradon(squeeze(block(:,k,1:end-1)),angles,'spline','Hann',1);
    
    end
    
% 
%     % generate variable name
%     v=strcat(nam_block,nam_stp{i},'_',pfx_cas{j});
%     eval([v '= blocklog;']);
%     %save block to disk
%     if ind_write==1
%         fprintf('saving %s to harddrive.....',v)
%         savefast(strcat(dir_work,v),v);
%         fprintf('done\n')
%     end
%     eval(['clear ' v]); % delete from workspace
%     end
% 
%     clear block blocklog cl3
%     t2=toc;
%     fprintf(t_string,name,t2-t1,t2/60);
end