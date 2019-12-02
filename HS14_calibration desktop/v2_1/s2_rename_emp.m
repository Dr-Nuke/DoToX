%process convention
% 1 - data
% 2 - 180°
% 3 - Open beam

%% empty
cas=1;
stp=1;
dir_r=strcat(dir_orig,dir_raw{cas});
dir_w=strcat(dir_work,dir_stp{stp},dir_cas{cas});


emp_xls_name  = 'empty.xls';
emp_xls_path  = strcat(dir_r,emp_xls_name);
emp_xls=f2_xls(emp_xls_path);



f2_boDir(dir_w)       
        
for i=nfiles{1}, %number of files %empty: 406
    disp(i)
    if mod(i,100)==0 %command line output
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
    elseif  
        fnn=sprintf('%04d',i); %file name number
        fnnn=sprintf('%04d',i); %file name number new
        name_fits=strcat('empty_',fnn,'.fits');
        name_fits_new=strcat('emp_',fnnn,'.fits');

        fileread=strcat(dir_r,name_fits);
        filewrite=strcat(dir_w,name_fits_new);

                
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

;%10
