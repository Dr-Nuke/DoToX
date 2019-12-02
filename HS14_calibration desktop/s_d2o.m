% import and sort the d2o channel measurement files


%% import the logfile data (.xls)
path_d2o      = strcat(readpath,'d2o\');
d2o_xls_name  = 'D2O.xls';
d2o_xls_path  = strcat(path_d2o,d2o_xls_name);
%empty_xls       = readtable(empty_xls_path);
load('d2oxls.mat');

for i=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
    d2o_log(:,i)=d2o_xls.(i);
end

%preallocate
d2o = zeros(length(xpix),length(ypix),max(empty_xls.(2))); %image data
d2o_OB = zeros(length(xpix),length(ypix)); %oenbeam
d2o_180= zeros(length(xpix),length(ypix)); %180° measurement


fprintf('\n D2O images skipped:');
%% read in full d2ochannel fits


names={'D2O_','D2O_r_'};
imax=[936,681];
idx_im=0; %image counter index

for j=1:2 % for the fact that when the experiment was continued after 
            % the break at i=936, i started over with 1
    
    for i=1:imax(j), 
        idx_im=idx_im+1;
        d2o_log(idx_im,13)=i;
        d2o_log(idx_im,14)=idx_im;
        
        % check if the image is correct
        ii=d2o_log(idx_im,2)+1; % projection name number
        fnn=sprintf('%04d',i); %file name number
        name_fits=strcat(names{j},fnn,'.fits');
        filepath=strcat(path_d2o,name_fits);

%         if mod(i,100)==0, % command line output
%             disp(sprintf('reading image %s', name_fits));
%         end
       
        if  any([d2o_log(idx_im,7)==0, d2o_log(idx_im,2)==375 ]),
            % = "image was bad" or "the 360° = 0° image"
            fprintf(' %d', i);
            d2o_log(idx_im,11:14)=nan;
            continue % bad image: ignore and goto next!
            
        elseif d2o_log(idx_im,1)==1 % actual d2ochannel image
       d2o(:,:,ii)=d2o(:,:,ii)+f_fits(filepath,thresh_medi,xpix,ypix);

        elseif d2o_log(idx_im,1)==3, %openbeam image
            d2o_OB =d2o_OB+ f_fits(filepath,thresh_medi,xpix,ypix);
            
        elseif d2o_log(idx_im,1)==2, %180° image
            d2o_180 =d2o_180+ f_fits(filepath,thresh_medi,xpix,ypix);
        end

        
        d2o_log(idx_im,11)=ii; % make a note in the log

        %timestamp 
        tmp=datevec(datenum(d2o_xls.(5){idx_im}))-datevec(datenum(d2o_xls.(4){idx_im}));
        d2o_log(idx_im,12)=mod(tmp(6)+60,60);
        d2o_log(idx_im,13)=i;
        


    end
end    
%% image corrections
% now sort out the bad images from the list and keep relevant data
d2o_log2=d2o_log([isfinite(d2o_log(:,12))],[2,6,10:14]);

%now deterimine the number of images per projection
[counts,centers] = hist(d2o_log2(1:1498,1),[0:374]);

%and divide each projection by the number of summed up images
d2o=d2o./repmat(reshape(counts(:),1,1,length(counts(:))),length(xpix),length(ypix));

d2o_OB=d2o_OB/4; % since its 4 images a 20s taken
d2o_180=d2o_180/4;
