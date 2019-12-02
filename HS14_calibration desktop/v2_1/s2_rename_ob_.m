%% ob

ob_file='ob_000';




ii=[1,3,4,5,6]; % see log file C:\data\tomo_HS14\01_apcs\New folder\ob.xls
for i =1:length(ii);  
    name_fits=strcat(ob_file,num2str(ii(i)),'.fits');
    name_fits_new=strcat('ob_',num2str(i),'.fits');

    fileread=strcat(readpath,'openbeam\',name_fits);
    filewrite=strcat(writepath,'1raw\ob_\',name_fits_new);

        if ind_write==1
            copyfile(fileread,filewrite)
        end
end






