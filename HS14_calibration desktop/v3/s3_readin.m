%% script to read the files into the matlab block

if rin
    if d2o
        
        % organize the log table
        dir_r=strcat(dir_orig,'d2o\');                  % read path
        f3_FolderCheck(dir_work)              %check and create write path       
        xls_path  = strcat(dir_r,'D2O.xls');
        d2o_tab_rin=f2_xls(xls_path);
        d2o_tab_rin.Properties.VariableNames{12}='slice_Nr';
        name_fits_new=' '; % parfor warning fix
        
        blocksize=850;
        ind_min=1:blocksize:1650;
        d2o_dat_rin=zeros(1200,2450,blocksize,'single');%sum(xls_tab{:,7})); % only as many as good images
        
        kk=1; %slice number counter
        
        for k = nfiles{2}{1}
            if any([d2o_tab_rin(k,7)==0, d2o_tab_rin(k,10)==360 ])
                % check if bad image or 360° redundant image => discard
                d2o_tab_rin(k,:)=[];
            end
        end
            
            
            
            

    end
end





%1

%2

%3

%5

%6

%7

%8

%9

%10