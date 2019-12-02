%% renaming script

%check if we shall run this script
if any(itr_stp==1)
    
    tic                                         %count time
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...'])             %command line output            
    
    i=1; % specifies step index of this file, is one
    
    
    for j = itr_cas_emp   % iterate emp, d2o, cl3, dc
        dir_r=strcat(dir_orig,dir_raw{j});                  % read path
        dir_w=strcat(dir_work,dir_stp{i},dir_cas{j});   % write path
        f3_FolderCheck(dir_w)              %check and create write path

        
        % xls logfile
        xls_path  = strcat(dir_r,nam_xls{j});
        xls_tab=f2_xls(xls_path);
        name_fits_new=' '; % parfor warning fix
        
        %parfor- specific fixes
        nfiles21end=nfiles{2}{1}(end); 
        pfx_casj=pfx_cas{j};
        pfx_cas5=pfx_cas{5};
        pfx_cas6=pfx_cas{6};
        
        
        for jj= 1:ind_jj(j); %special D2O double iteration
            
            nam_raw2=nam_raw{j}{jj}; %parfor-specific fix
            
            for k=nfiles{j}{jj}
                f3_BoCountDisp(k,300);  % command line counter
                xls_work=xls_tab(k,:);  % slice required for parfor

                if xls_work{1,1}==1 % if it is a data projection 
                    xls_work{1,2}=xls_work{1,2}+1; %increase the projection number
                end                                %such that it starts at 1

                %check for validity of image
                if  any([xls_work{1,7}==0, xls_work{1,2}==376 ]), %if not valid, then...
                    % = "image was bad" or "the 360° = 0° image"
                    % then skip
                    
%                     if xls_work{1,7}==1 %set the 'good' property
%                         xls_work{1,7}=0;% to 'bad' if not yet so
%                     end

                elseif  xls_work{1,7}==1
                    if and(j==1,xls_work{1,1}==2) % ob images for empty
                        xls_work{1,1}=3;
                    end
                    
                    %raw file name
                    fnn=sprintf('%04d',k); %file name number
                    if jj==1 %d20 first case
                        name_fits=strcat(nam_raw2,fnn,'.fits'); %creae r-file name string
                    else    %d20 special case
                        fnnn=sprintf('%04d',k-nfiles21end); %file name number
                        name_fits=strcat(nam_raw2,fnnn,'.fits'); %creae r-file name string
                    end
                    
                    
                    % tabel entries and w-file names
                    if xls_work{1,1}==1 %data image
                        name_fits_new=strcat(pfx_casj,fnn,'_','.fits');
                        xls_work{1,12}={'data'};
                    elseif xls_work{1,1}==2 % 180+ image
                        name_fits_new=strcat(pfx_casj,fnn,'_',pfx_cas6,'.fits');
                        xls_work{1,12}={pfx_cas6};
                    elseif xls_work{1,1}==3 % OB° image
                        name_fits_new=strcat(pfx_casj,fnn,'_',pfx_cas5,'.fits');
                        xls_work{1,12}={pfx_cas6};
                    else % wtf?!?
                        disp([i,j,k,jj,kk]);
                        error('image category number is bad');
                        
                        %parfor-specific warning message suppression 
                        % (try & outcomment this one line & see what
                        % happens)                       
                        name_fits_new=' ';
                    end
                    fileread=strcat(dir_r,name_fits);
                    filewrite=strcat(dir_w,name_fits_new);
                    xls_work{1,13}={name_fits_new};
                    
                    %finally copy file
                    if ind_write
                        copyfile(fileread,filewrite)
                    end
                end %if
                
                xls_tab(k,:)=xls_work;
                
            end % parfor
            
        end % special d2o
        v = genvarname(['xls_',pfx_cas{j}]);
        eval([v '=xls_tab;']); %wite down  table to variable
        
    end % j-loop
    
    toc 
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

