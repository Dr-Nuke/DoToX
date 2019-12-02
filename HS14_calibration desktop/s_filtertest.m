
%close all
path_d2o      = strcat(path,'d2o\');
ii=[674:677];
thresh=0.02;
m=3; %figure number
n=4;

n_im=[1,2;3,4;5,6]
medkern=[3,3;1,3;3,1]

for i=1:3,
    %reand in image
    
    fnn=676; %file name number
    path=strcat(readpath,'D2O\D2O_r_0676.fits');
    tmp=fitsread(path)';
    im=tmp(250:650,1:400);
    im_(:,:,i)=im;
    
    %plot
    figure(n_im(i,1))
    title(ii(i))
    subplot(2,2,1)
    imbo4(im); caxis([0,55000])
    
    %apply filter
    filt=medfilt2(im,medkern(1,:));
    filt_(:,:,i)=filt;
    diff=im-filt;
    diff_(:,:,i)=diff; 
    
    
    %calculate correctin
    corr=im;
    corr(diff>thresh*im)=filt(diff>thresh*im);
    corr_(:,:,i)=corr;
    subplot(2,2,3); imbo4(corr.*(diff>thresh*im)); caxis([0,55000])
    
    
    
    %plot correction
    subplot(2,2,4)
    imbo4(f_medifilter1(im,thresh)); caxis([0,55000])
    
    
     truesize;



%%
figure(n_im(i,2))
threshh=[0.01,0.02,0.03,0.04,0.05,0.06];

for j =1:6
    subplot(2,3,j)
    imbo4(corr.*(diff>threshh(j)*im)); caxis([0,64000]);
    [a,b]=size(corr);
    n=squeeze(sum(sum(diff>threshh(j)*im)));
    titlestring=sprintf('thresh= %f, frac= %f',threshh(j),n/(a*b));
   title(titlestring);
  nn(i+1,j)=n/(a*b);
end

truesize()
end
nn(1,:)=threshh;
disp(nn)
    

%%
nbin=1000;
figure(8)
hist(im(:),nbin);
figure(9)
a=f_medifilter1(im,0.02);
hist(a(:),nbin);
figure(10)
b=im-a;
hist(b(:),nbin); %ylim([0,10])
    
    
    
    
    
    
    
    
    
    
    
    