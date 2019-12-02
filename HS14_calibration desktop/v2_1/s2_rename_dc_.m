%% DC

dc_file='Boratedpoly_block_';

imax=5;

for i =1:imax,
    name_fits=strcat(dc_file,num2str(i),'.fits');
    name_fits_new=strcat('dc_',num2str(i),'.fits');

    fileread=strcat(readpath,'CD\',name_fits);
    filewrite=strcat(writepath,'1raw\dc_\',name_fits_new);

        if ind_write==1
            copyfile(fileread,filewrite)
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