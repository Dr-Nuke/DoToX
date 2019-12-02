%% this script reconstructs the image

n_rec=749; % number of pixels the reconstuction is high & wide
path_cl3_rec=strcat(writepath,'CL3\');
path_d2o_rec=strcat(writepath,'d2o\');
path_emp_rec=strcat(writepath,'empty\');

%mask
[x,y]=meshgrid(1:n_rec,1:n_rec);
mask=zeros(n_rec);
mask((x-floor(0.5*n_rec)).^2+(y-floor(0.5*n_rec)).^2<floor(0.49*n_rec)^2)=1;

% centering for each plane
%for i=[1:length(ypix)] %1:length(ypix)
for i=[50:60]    
    
    % do the cross correlations
    if ind_debug==1
        [cl3_c,lags_cl3]=xcov(squeeze(raw_cl3_d2o(:,i,1)),flipud(raw_cl3_d2o_180(:,i)));
        [d2o_c,lags_d2o]=xcov(squeeze(raw_d2o_emp(:,i,1)),flipud(raw_d2o_emp_180(:,i)));
        [emp_c,lags_emp]=xcov(squeeze(empty_corr(:,i,1)),flipud(empty_180_corr(:,i)));
    elseif ind_debug==0
        [cl3_c,lags_cl3]=xcov(squeeze(cl3(:,i,1)),flipud(cl3_180(:,i)));
        [d2o_c,lags_d2o]=xcov(squeeze(d2o(:,i,1)),flipud(d2o_180(:,i)));
        [emp_c,lags_emp]=xcov(squeeze(empty(:,i,1)),flipud(empty_180(:,i)));
    else
        disp('bad ind_debug');
    end
        
    %
    %    % here fill in sub-pixel adjustment
    %
    
    
    % doing the max lag and its shift value
    [cl3_maxlag(i,1),cl3_maxlag(i,2)]=max(cl3_c);
    [d2o_maxlag(i,1),d2o_maxlag(i,2)]=max(d2o_c);
    [emp_maxlag(i,1),emp_maxlag(i,2)]=max(emp_c);
 
    cl3_maxlag(i,3)=lags_cl3(cl3_maxlag(i,2));
    d2o_maxlag(i,3)=lags_d2o(d2o_maxlag(i,2));
    emp_maxlag(i,3)=lags_emp(emp_maxlag(i,2));
    
    % actually centering
    if ind_debug==1
        cl3_cent{i}=f_centCrop(squeeze(squeeze(raw_cl3_d2o(:,i,:))),cl3_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
        d2o_cent{i}=f_centCrop(squeeze(squeeze(raw_d2o_emp(:,i,:))),d2o_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
        emp_cent{i}=f_centCrop(squeeze(squeeze(empty_corr(:,i,:))),emp_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
    elseif ind_debug==0
        cl3_cent{i}=f_centCrop(squeeze(squeeze(cl3(:,i,:))),cl3_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
        d2o_cent{i}=f_centCrop(squeeze(squeeze(d2o(:,i,:))),d2o_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
        emp_cent{i}=f_centCrop(squeeze(squeeze(empty(:,i,:))),emp_maxlag(i,3),xmin,xmax,xdo_1,xdo_2);
    else
        disp('bad ind_debug');
    end
    


  
    % reconstruction including rotation by 12 steps such that the cooling
    % channels are foind in the NE, NW, SW ad SE corner
    d2o_rec(:,:,i)=iradon(circshift(d2o_cent{i},12,2),angles,'spline','Han',1,n_rec);%.*mask;
    cl3_rec(:,:,i)=iradon(circshift(cl3_cent{i},12,2),angles,'spline','Han',1,n_rec);%.*mask;
    emp_rec(:,:,i)=iradon(circshift(emp_cent{i},12,2),angles,'spline','Han',1,n_rec);%.*mask;
    
    fitswrite(d2o_rec(:,:,i),strcat(path_d2o_rec,sprintf('D2O_%04d',ymin+i-1),'.fits'))
    fitswrite(cl3_rec(:,:,i),strcat(path_cl3_rec,sprintf('CL3_%04d',ymin+i-1),'.fits'))
    fitswrite(emp_rec(:,:,i),strcat(path_emp_rec,sprintf('Empty_%04d',ymin+i-1),'.fits'))
    
    disp([i,cl3_maxlag(i,3),d2o_maxlag(i,3),emp_maxlag(i,3)])
    
%     cl3_rec2=iradon(-log(cl3_test-d2o_test),angles,'spline','Han',1,n_rec);
%     cl3_rec3=iradon(-log(cl3_test-empty_test),angles,'spline','Han',1,n_rec);
%     d2o_rec2=iradon(-log(d2o_test-empty_test),angles,'spline','Han',1,n_rec);
% 
%     fitswrite(cl3_rec2,strcat(path_cl3_rec,sprintf('2_%04d',i),'.fits'))
%     fitswrite(cl3_rec3,strcat(path_cl3_rec,sprintf('3_%04d',i),'.fits'))
%     fitswrite(d2o_rec2,strcat(path_d2o_rec,sprintf('2_%04d',i),'.fits'))


    
    
end

