% import and sort the cl3 channel measurement files


%% import the logfile data (.xls)
path_cl3      = strcat(readpath,'chcl3\');
cl3_xls_name  = 'CHCL3.xls';
cl3_xls_path  = strcat(path_cl3,cl3_xls_name);
%empty_xls       = readtable(empty_xls_path);
load('cl3xls.mat');

for i=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
    cl3_log(:,i)=cl3_xls.(i);
end




%% read in full cl3channel fits
%preallocate

cl3 = zeros(length(xpix),length(ypix),max(empty_xls.(2))); %image data
cl3_OB = zeros(length(xpix),length(ypix)); %oenbeam
cl3_180= zeros(length(xpix),length(ypix)); %180° measurement


idx_im=0; %image counter index

fprintf('\n Cl3 images skipped:');

for i=1:1621
    idx_im=idx_im+1;
    cl3_log(idx_im,14)=i;
    cl3_log(idx_im,14)=idx_im;

    % check if the image is correct
    ii=cl3_log(idx_im,2)+1; % projection name number
    fnn=sprintf('%04d',i); %file name number
    name_fits=strcat('CHCL3_',fnn,'.fits');
    filepath=strcat(path_cl3,name_fits);
    
%     if mod(i,100)==0, % command line output
%         disp(sprintf('reading image %s', name_fits));
%     end

    if  any([cl3_log(idx_im,7)==0, cl3_log(idx_im,2)==375 ]),
        % = "image was bad" or "the 360° = 0° image"
        fprintf(' %d', i);
        cl3_log(idx_im,11:14)=nan;
        continue % bad image: ignore and goto next!
    elseif cl3_log(idx_im,1)==1 % actual cl3channel image
        cl3(:,:,ii)=cl3(:,:,ii)+f_fits(filepath,thresh_medi,xpix,ypix);

    elseif cl3_log(idx_im,1)==3, %openbeam image
        cl3_OB =cl3_OB+ f_fits(filepath,thresh_medi,xpix,ypix);

    elseif cl3_log(idx_im,1)==2, %180° image
        cl3_180 =cl3_180+ f_fits(filepath,thresh_medi,xpix,ypix);  
    end

    cl3_log(idx_im,11)=ii; % make a note in the log

    %timestamp 
    tmp=datevec(datenum(cl3_xls.(5){idx_im}))-datevec(datenum(cl3_xls.(4){idx_im}));
    cl3_log(idx_im,12)=mod(tmp(6)+60,60);
    cl3_log(idx_im,13)=i;



end
  
%% 
% now sort out the bad images from the list and keep relevant data
cl3_log2=cl3_log([isfinite(cl3_log(:,12))],[2,6,10:14]);

%now deterimine the number of images per projection
[counts,centers] = hist(cl3_log2(1:1500,1),[0:374]);

%and divide each projection by the number of summed up images
cl3=cl3./repmat(reshape(counts(:),1,1,length(counts(:))),length(xpix),length(ypix));

cl3_OB=cl3_OB/4; % since its 4 images a 20s taken
cl3_180=cl3_180/4;
