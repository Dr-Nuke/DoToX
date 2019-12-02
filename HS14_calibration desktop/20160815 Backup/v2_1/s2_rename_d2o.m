%process convention
% 1 - data
% 2 - 180°
% 3 - Open beam%% D20

%% d2o

path_d2o      = strcat(readpath,'D2O\');
d2o_xls_name  = 'D2O.xls';
d2o_xls_path  = strcat(path_d2o,d2o_xls_name);
d2o_xls=f2_xls(d2o_xls_path);


names={'D2O_','D2O_r_'};
imax=[938,681]; % ATTENTION! I added files 937 and 938 as they are listed
% in the .xls but not existent.



idx_im=0; %image counter index

for j=1:2 % for the fact that when the experiment was continued after 
            % the break at i=936, i started over with 1

    for i=1:imax(j),
        idx_im=idx_im+1;
        %check if the image is correct  

        if mod(idx_im,100)==0
             fprintf('_%d %d',i,idx_im) 
        end

        if d2o_xls{idx_im,1}==1 %increase the projection number
            d2o_xls{idx_im,2}=d2o_xls{idx_im,2}+1;    %such that it starts at 1
        end
        
        % write file numbers
        d2o_xls{idx_im,13}=i;
        d2o_xls{idx_im,14}=idx_im;

        %check for validity of image
        if  any([d2o_xls{idx_im,7}==0, d2o_xls{idx_im,2}==376 ]),
            % = "image was bad" or "the 360° = 0° image"
            % then skip
            if d2o_xls{idx_im,7}==1 %ste the 'good' property
                d2o_xls{idx_im,7}=0;
            end
            continue
        else 
            fnn=sprintf('%04d',i); %file name number
            fnnn=sprintf('%04d',idx_im); %file name number new
            name_fits=strcat(names{j},fnn,'.fits');
            name_fits_new=strcat('d2o_',fnnn,'.fits');
            
            fileread=strcat(readpath,'D2O\',name_fits);
            filewrite=strcat(writepath,'1raw\d2o\',name_fits_new);
        
        


            if ind_write==1
                copyfile(fileread,filewrite)
            end
                
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