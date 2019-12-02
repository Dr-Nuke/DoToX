%% corrects all images according to xls file
load('xls_file')
flag={'emp';'d2o';'cl3'};
load('RecProp_file.mat');

dcfile=strcat('dc','_add_','.mat');
dcpath=strcat(writepath,'3add\','dc_','\');
doseXrange=[1:250,950:1120];
dc=f2_boLoad(strcat(dcpath,dcfile));

cropx=[231:970];

RecProp.corr_image_x=length(cropx);
RecProp.corr_image_y=RecProp.filtered_image_y;

save('RecProp_file.mat','RecProp');

for i =1:3 %number of falgs
    disp(strcat('correcting_ ',flag{i}))
    rpath=strcat(writepath,'3add\',flag{i},'\');
    wpath=strcat(writepath,'4cor\',flag{i},'\');

    
    for j=1:375 %number of angles of useful images
        disp(j)
        rfile=strcat(flag{i},'_add_',sprintf('%04d',j),'.mat');
        wfile=strcat(flag{i},'_cor_',sprintf('%04d',j),'.mat');
        
        im=f2_boLoad(strcat(rpath,rfile));
        im=im-dc; %dc correction
        im=im/f2_doseFactor(im,doseXrange); %dose normalization
        
        if ind_write==1
            f2_parforSave(strcat(wpath,wfile),im(cropx,:))
        end
        
        
    end
    
    % 180° images
    rfile=strcat(flag{i},'_add_','180_','.mat');
    wfile=strcat(flag{i},'_cor_','180_','.mat');
    
    im=f2_boLoad(strcat(rpath,rfile));
    im=im-dc; %dc correction
    im=im/f2_doseFactor(im,doseXrange); %dose normalization
    if ind_write==1
        f2_parforSave(strcat(wpath,wfile),im(cropx,:))    
    end
    
    % OB files
    rfile=strcat(flag{i},'_add_','OB__','.mat');
    wfile=strcat(flag{i},'_cor_','OB__','.mat');
    
    im=f2_boLoad(strcat(rpath,rfile));
    im=im-dc; %dc correction
    im=im/f2_doseFactor(im,doseXrange); %dose normalization
    if ind_write==1
        f2_parforSave(strcat(wpath,wfile),im(cropx,:))   
    end
end




        


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

