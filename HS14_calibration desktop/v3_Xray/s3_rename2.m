%% renaming script

%check if we shall run this script
if  any(itr_stp==1) % execute this script?
    t1=toc;
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...'])             %command line output            
    
    i=1; % specifies step index of this file, is one
    
    % set up logfile
    if exist(nam_log_ful,'file')==2
        load(nam_log_ful); %load logfile.mat, variable name is logfile
    else
        logfile={};
    end
    
    for j = itr_cas   % iterate emp, d2o, cl3, dc
        
        if j==4 %dark current extrawurst
            %% dark current
            dir_r=strcat(dir_orig,'CD/');      % read path
            dir_w=strcat(dir_work,dir_stp{i},'dc_/');   % write path
            f3_FolderCheck(dir_w)
            
            %make xls
            xls_path  = strcat(dir_r,'DC_fake.xls');
            xls_tab=f3_xls(xls_path);            

            for k=1:5
                file_r=strcat('Boratedpoly_block_',num2str(k),'.fits');
                fileread=strcat(dir_r,file_r);
                file_w=strcat('DC_',num2str(k),'.fits');
                filewrite=strcat(dir_w,file_w);
                if ind_write==1
                    copyfile(fileread,filewrite)
                end
                
                xls_tab{k,9}=k;     %slice number
                xls_tab{k,13}={file_r};     %old file name
                xls_tab{k,14}={file_w};     %new  file name


                
            end
        else
            
     
            dir_r=strcat(dir_orig,dir_raw{j});                  % read path
            dir_w=strcat(dir_work,dir_stp{i},dir_cas{j});   % write path
            f3_FolderCheck(dir_w)              %check and create write path


            % xls logfile
            % xls first column: 1 = data, 2= 180�, 3=ob; original empty 2=ob
            xls_path  = strcat(dir_r,nam_xls{j});
            xls_tab=f3_xls(xls_path);
            

            %parfor- specific fixes
            nfiles21end=nfiles{2}{1}(end); 
            pfx_casj=pfx_cas{j};
%             pfx_cas5=pfx_cas{5}; %removed due to error
%             pfx_cas6=pfx_cas{6};

            k=1; % slice number
            k_180=1; % 180� image number

            for jj= 1:ind_jj(j) %special D2O double iteration


                for kk=nfiles{j}{jj} %file number (not file name number!)
                    f3_BoCountDisp(k,300);  % command line counter

                    if xls_tab{k,1}==1 % if it is a data projection 
                        xls_tab{k,2}=xls_tab{k,2}+1; %increase the projection number
                    end                                %such that it starts at 1

                    %check for non-validity of image
                    if  any([xls_tab{k,7}==0,...    % image is faulty
                            xls_tab{k,2}==376,...   % redundant 360� images 
                            xls_tab{k,1}==3,...     % OB images for d2o and cl3
                            and(j==1,...
                                xls_tab{k,1}==2)])  % OB images for empty

                        % then delete the row and discard
                        xls_tab(k,:)=[];
                        % and dont increment k

                    elseif  xls_tab{k,7}==1         % good image

                        %raw file name
                        fnn=sprintf('%04d',kk); %old file name number
                        fnn_new=sprintf('%04d',k); %new file name number

                        if jj==1 %d20 first case
                            name_fits=strcat(nam_raw{j}{jj},fnn,'.fits'); %creae r-file name string
                        else    %d20 special case
                            fnnn=sprintf('%04d',kk-nfiles21end); %file name number
                            name_fits=strcat(nam_raw{j}{jj},fnnn,'.fits'); %creae r-file name string
                        end



                        % tabel entries and w-file names
                        if xls_tab{k,1}==1 %data image
                            name_fits_new=strcat(pfx_casj,fnn_new,'_','.fits');
                            xls_tab{k,12}={'data'};

                        elseif xls_tab{k,1}==2 % 180� image
                            name_fits_new=strcat(pfx_casj,'_180',num2str(k_180),'_.fits');
                            xls_tab{k,12}={'180'};
                            k_180=k_180+1;
                        else % wtf?!?
                            disp([i,j,k,jj,kk]);
                            error('image category number is bad');

                            %parfor-specific warning message suppression 
                            % (try & outcomment this one line & see what
                            % happens)                       
                            % name_fits_new=' ';
                        end

                        xls_tab{k,9}=k;     %slice number
                        xls_tab{k,13}={name_fits};     %old file name
                        xls_tab{k,14}={name_fits_new};     %new  file name

                        fileread=strcat(dir_r,name_fits);
                        filewrite=strcat(dir_w,name_fits_new);


                        %finally copy file
                        if ind_write==1
                            copyfile(fileread,filewrite)
                        end
                        k=k+1;
                    end %if



                end % parfor

            end % special d2o
        end
        
        if j==1 % empty 180 deg manual fix
            
            
            
            file_r='empty_000_180.fits';
            fileread=strcat(dir_r,file_r);
            
            file_w='emp__1801.fits';
            filewrite=strcat(dir_w,file_w);
            if ind_write==1
                copyfile(fileread,filewrite)
            end
            xls_tab(end+1,:)=xls_tab(end,:);
            xls_tab(end,1)={2};
            xls_tab(end,2)={0};
            xls_tab(end,9)={376};
            xls_tab(end,10)={180};
            xls_tab{end,12}={'180_'};
            xls_tab{end,13}={file_r};
            xls_tab{end,14}={file_w};            
        end
            
     % write down log file
     logfile{i}{j}=xls_tab;
      
    end % j-loop
    %save logfile to disk
    save(nam_log_ful,'logfile');
    
    t2=toc
    fprintf(t_string,name,t2-t1,t2/60)
end %if





 %1

 %2

 %3

%5

%6

%7

%8



%9

%10



















%

