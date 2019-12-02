%% cl3
%process convention
% 1 - data
% 2 - 180°
% 3 - Open beam

%% CL3

path_cl3      = strcat(readpath,'chcl3\');
cl3_xls_name  = 'chcl3.xls';
cl3_xls_path  = strcat(path_cl3,cl3_xls_name);
cl3_xls=f2_xls(cl3_xls_path);


parfor i=1:nfiles(3), %number of files %cl3 1621
    
    % command window output
    if mod(i,100)==0
         fprintf('_%d ',i) 
    end

    if cl3_xls{i,1}==1 %increase the projection number
        cl3_xls{i,2}=cl3_xls{i,2}+1;    %such that it starts at 1
    end
    
    % write file numbers
    cl3_xls{i,13}=i;
    cl3_xls{i,14}=i;
    
    %check for validity of image        
    if  any([cl3_xls{i,7}==0, cl3_xls{i,2}==376 ]),
        % = "image was bad" or "the 360° = 0° image"
        % then skip
        if cl3_xls{i,7}==1 %ste the 'good' property
            cl3_xls{i,7}=0;
        end
        continue
    else 
        fnn=sprintf('%04d',i); %file name number
        fnnn=sprintf('%04d',i); %file name number new
        name_fits=strcat('CHCL3_',fnn,'.fits');
        name_fits_new=strcat('cl3_',fnnn,'.fits');

        fileread=strcat(readpath,'chcl3\',name_fits);
        filewrite=strcat(writepath,'1raw\cl3\',name_fits_new);
        

        if ind_write==1
            copyfile(fileread,filewrite)
        end
        
    end
end

fprintf('\n') 

; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

; %10