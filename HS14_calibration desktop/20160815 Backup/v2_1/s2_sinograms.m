%% make the sinograms

%% corrects all images according to xls file
load('xls_file')
flag={'emp';'d2o';'cl3'};
load('RecProp_file.mat');

x=RecProp.corr_image_x;
y=RecProp.corr_image_y;
ang=RecProp.numAngles;


for i =1:3 %number of flags
    disp(strcat('creating_sin_ ',flag{i}))
    rpath=strcat(writepath,'4cor\',flag{i},'\');
    wpath=strcat(writepath,'5sin\',flag{i},'\');
    
    f2_FolderCheck({wpath,rpath})
    
    BigData=zeros(ang,x,y); 
    
    for j=1:ang %number of angles of useful images

        rfile=strcat(flag{i},'_cor_',sprintf('%04d',j),'.mat');
        BigData(j,:,:)=f2_boLoad(strcat(rpath,rfile));
    end
        
    disp('writing sinogram for plane ')          
    for j=1:y
        if mod(j,200)==0
            fprintf('%i',j)
        end
        wfile=strcat(flag{i},'_sin_',sprintf('%04d',j),'.mat');
        sin=squeeze(BigData(:,:,j));
        if ind_write==1
            f2_parforSave(strcat(wpath,wfile),sin) 
        end
    end
    
        
    save(strcat(wpath,strcat('TotalSino_',flag{i},'.mat')),'BigData','-v7.3') 
    
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


