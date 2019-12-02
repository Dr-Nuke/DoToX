% import and sort the empty channel measurement files


%% import the logfile data (.xls)
path_empty      = strcat(readpath,'empty\');
empty_xls_name  = 'empty.xls';
empty_xls_path  = strcat(path_empty,empty_xls_name);
%empty_xls       = readtable(empty_xls_path);
load('empxls.mat');

for i=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
    empty_log(:,i)=empty_xls.(i);
end


fprintf('\n Empty images skipped:');

%% read in full emptychannel fits
%preallocate
empty = zeros(length(xpix),length(ypix),max(empty_xls.(2)));
empty_OB = zeros(length(xpix),length(ypix));

for i=1:406, %number of files %empty: 406
    % check if the image is correct
    ii=empty_log(i,2)+1; % "correct projection" name number
    fnn=sprintf('%04d',i); %file name number
    name_fits=strcat('empty_',fnn,'.fits');
    filepath=strcat(path_empty,name_fits);
    
%     if mod(i,50)==0, % command line output
%         disp(sprintf('reading image %s', name_fits));
%     end
    
    if  any([empty_log(i,7)==0, empty_log(i,2)==375 ]) ,
        % =or(image was bad or the 360° = 0° image it is
        fprintf('%d ', i);
        empty_log(i,11:12)=nan;
        continue % bad image: ignore and goto next!
        
    elseif empty_log(i,1)==2, %openbeam image
        
        empty_OB = empty_OB+ f_fits(filepath,thresh_medi,xpix,ypix);
    else % actual emptychannel image
        empty(:,:,ii)=empty(:,:,ii)+f_fits(filepath,thresh_medi,xpix,ypix);
    end
    
    empty_log(i,11)=ii; % make a note in the log
    
    %timestamp 
    tmp=datevec(datenum(empty_xls.(5){i}))-datevec(datenum(empty_xls.(4){i}));
    empty_log(i,12)=mod(tmp(6)+60,60);
    empty_log(i,13)=i;
    
    
end
%%
empty_OB=empty_OB/5; % since its 5 images a 20s taken for empty
    
% now sort out the bad images from the list and keep relevant data
empty_log2=empty_log([isfinite(empty_log(:,12))],[2,6,10,11,12,13]);

%manually read in 180° image
empty_180=f_fits(strcat(path_empty,'empty_000_180.fits'),thresh_medi,xpix,ypix);


