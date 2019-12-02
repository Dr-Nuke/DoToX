%% pre-filtering script
%  DC preparing






%% actual filtering
%check if we shall run this script
if  any(itr_stp==3) % execute this script?
    i=3; % specifies step index of this file, is two
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
    load(nam_log_ful) %load logfile.mat, variable name is logfile
    j=4
    log_wrk=logfile{i-1}{j}; %working copy of the logfile
    n_blocks=max(log_wrk{:,7}); % number of blocks for this case

    % k is the image number
    % kk is the block number
    % kkk i the index number of image k within block kk 

    for kk=1:n_blocks
        disp(strcat(nam_block,nam_stp{i},'_',pfx_cas{j},num2str(kk)));
        %load the block from previous step
        block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{j},num2str(kk));
        path_block=strcat(dir_work,block_old);
        DC=f3_loadSingleVariableMATFile(path_block);
        [x,y,z]=size(DC);
        
        %find the amount of images in that block



        for k=1:z
            
            im=DC(:,:,k);
            %filter for bright & dark spots
            im=f3_MediFilter2(im,kernel_medi2,thresh_medi2,'median' );

            % filter gammaspots(median)
            % after-filtering to remove circles left overy by spotfilter
            im=f_MediFilter1(im,thresh_medi1);    %vertical
            im=f_MediFilter1_2(im,thresh_medi1);  % horizontal (or vice versa)  

            DC(:,:,k)=im;


        end %k-loop ends
    end
    
    DCpath=strcat(dir_work,'DC');
    save(DCpath,'DC');
    
    

    
end