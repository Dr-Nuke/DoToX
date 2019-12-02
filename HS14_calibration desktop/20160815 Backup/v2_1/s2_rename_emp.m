%process convention
% 1 - data
% 2 - 180°
% 3 - Open beam

%% empty

path_emp      = strcat(readpath,'empty\');
emp_xls_name  = 'empty.xls';
emp_xls_path  = strcat(path_emp,emp_xls_name);
emp_xls=f2_xls(emp_xls_path);


for i=1:406, %number of files %empty: 406
    if mod(i,100)==0
         fprintf('_%d ',i) 
    end

    if emp_xls{i,1}==1 %increase the projection number
        emp_xls{i,2}=emp_xls{i,2}+1;    %such that it starts at 1
    end
    
    % write file numbers
    emp_xls{i,13}=i;
    emp_xls{i,14}=i;

    %check for validity of image
    if  any([emp_xls{i,7}==0, emp_xls{i,2}==376 ]),
        % = "image was bad" or "the 360° = 0° image"
        % then skip
        if emp_xls{i,7}==1 %ste the 'good' property
            emp_xls{i,7}=0;
        end
        continue
    else 
        fnn=sprintf('%04d',i); %file name number
        fnnn=sprintf('%04d',i); %file name number new
        name_fits=strcat('empty_',fnn,'.fits');
        name_fits_new=strcat('emp_',fnnn,'.fits');

        fileread=strcat(readpath,'empty\',name_fits);
        filewrite=strcat(writepath,'1raw\emp\',name_fits_new);

                
        if ind_write==1
            copyfile(fileread,filewrite)
        end
        
        if emp_xls{i,1}==2 % open beam: should be 3
            emp_xls{i,1}=3;   
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