%% adds all images according to xls file
load('xls_file')
flag={'emp';'d2o';'cl3'};

for i =1:3 %number of falgs
    disp(strcat('adding_ ',flag{i}))
    rpath=strcat(writepath,'2flt\',flag{i},'\');
    wpath=strcat(writepath,'3add\',flag{i},'\');
    
    parfor j=1:375 %number of angles of useful images
        
        subtable=xls{i}(xls{i}{:,2}==j.*xls{i}{:,7}==1,:);
        addIm=zeros(1120,2050);
        kmax=size(subtable,1);
        for k=1:kmax
            rfile=strcat(flag{i},'_filt_',sprintf('%04d',subtable{k,14}),'.mat');
            tmp=load(strcat(rpath,rfile));
            addIm=addIm+tmp.var/kmax;

        end
        wfile=strcat(flag{i},'_add_',sprintf('%04d',j),'.mat');
        if ind_write==1
            s2_parforSave(strcat(wpath,wfile),addIm)
        end
        
    end
    % open beams
    subtable=xls{i}(xls{i}{:,1}==3.*(xls{i}{:,7}==1),:);
    addIm=zeros(1120,2050);
    kmax=size(subtable,1);
    for k=1:kmax
        rfile=strcat(flag{i},'_filt_',sprintf('%04d',subtable{k,14}),'.mat');
        tmp=load(strcat(rpath,rfile));
        addIm=addIm+tmp.var/kmax;

    end
    wfile=strcat(flag{i},'_add_','OB__','.mat');
    if ind_write==1
        s2_parforSave(strcat(wpath,wfile),addIm)
    end
    
        % 180° images
    if any([2,3]==i)    % not for empty, they are extra
        % indices of images that are  180°,  show channel       and  are 'good'
        subtable=xls{i}(xls{i}{:,10}==180.*(xls{i}{:,8}==-19.99).*(xls{i}{:,7}==1),:);
        addIm=zeros(1120,2050);
        kmax=size(subtable,1);
        for k=1:kmax
            rfile=strcat(flag{i},'_filt_',sprintf('%04d',subtable{k,14}),'.mat');
            tmp=load(strcat(rpath,rfile));
            addIm=addIm+tmp.var/kmax;

        end
        wfile=strcat(flag{i},'_add_','180_','.mat');
        if ind_write==1
            s2_parforSave(strcat(wpath,wfile),addIm)
        end
        
    elseif i=1   % empty 180 is extra
        rpath=strcat(writepath,'2flt\',flag{i},'\');
        rfile=strcat(flag{1},'_filt_','180_','.mat');
        wfile=strcat(flag{i},'_add_','180_','.mat');
        
        im=im2double(fitsread(strcat(rpath,rfile))');
        if ind_write==1
            s2_parforSave(strcat(wpath,wfile),addIm)
        end
    end
    
end

%% DC images
        
disp(strcat('adding ','dc'))
addIm=zeros(1120,2050);
kmax=5;
rpath=strcat(writepath,'2flt\','dc_','\');
for k=1:kmax
    rfile=strcat('dc','_filt',num2str(k),'.mat');
    tmp=load(strcat(rpath,rfile));
    addIm=addIm+tmp.var/kmax;

end
wfile=strcat('dc','_add_','.mat');
wpath=strcat(writepath,'3add\','dc_','\');
if ind_write==1
    s2_parforSave(strcat(wpath,wfile),addIm)
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

; %10;
