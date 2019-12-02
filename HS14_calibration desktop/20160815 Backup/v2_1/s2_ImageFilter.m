
%% 1) emp
% cl3 special image errors start at 1104 (lines) and 1552 (wiered)
load('xls_file')
flag={'emp';'d2o';'cl3'};
for i =1:size(flag) %iterate cases
    disp(sprintf('filtering %s',flag{i}))
    jmax=size(xls{i},1);
    rpath=strcat(writepath,'1raw\',flag{i},'\');
    wpath=strcat(writepath,'2flt\',flag{i},'\');
    
    parfor j=1:jmax
        %check if the file is good
        if xls{i}{j,7}~=0

            rfile=strcat(flag{i},'_',sprintf('%04d',j),'.fits');
            wfile=strcat(flag{i},'_filt_',sprintf('%04d',j),'.mat');
            im=im2double(fitsread(strcat(rpath,rfile))');

            % crop away area of bad DC image
            % border pixel calculation

            im=im(81:end,201:2450);
            
            %line filter for special camera lines
            if and(i==3,j>=1104) %if cl3 only starting at 1104
                im=f_CamLineFilt(im);
            end

            %filter for bright & dark spots
            im=f_MediFilter2( im,kernel_medi2,thresh_medi2,'median' );

            % filter gammaspots(median)
            % after-filtering to remove circles left overy by spotfilter
            im=f_MediFilter1(im,thresh_medi1);
            im=f_MediFilter1_2(im,thresh_medi1);
            
            im=im(6:end-5,6:end-5); %cout out filter's bad area
            if ind_write==1
                s2_parforSave(strcat(wpath,wfile),im)  % parfor-proof save func
            end

            if mod(j,100)==0
                fprintf('%d ',j)
            end
            
        end
    end
end
    
%% DC

dc_file='Boratedpoly_block_';

imax=5;
rpath=strcat(writepath,'1raw\dc_\');
wpath=strcat(writepath,'2flt\dc_\');
disp('filtering DC')

parfor i =1:imax,
    rfile=strcat('dc_',num2str(i),'.fits');
    wfile=strcat('dc_filt',num2str(i),'.mat');
            im=im2double(fitsread(strcat(rpath,rfile))');

            % crop away area of bad DC image
            % border pixel calculation

            im=im(81:end,201:2250);

            % filter gammaspots(median)
            im=f_MediFilter1(im,thresh_medi1);

            %filter for bright & dark spots
            im=f_MediFilter2( im,kernel_medi2,thresh_medi2,'median' );
            
            im(6:end-5,6:end-5); %cout out filter's bad area
            if ind_write==1
                s2_parforSave(strcat(wpath,wfile),im)  % parfor-proof save func
            end
end

%% 180 empty
rpath=('C:\data\tomo_HS14\02_rawdata\empty\');
rfile=strcat('empty_000_180.fits');
wfile=strcat(flag{1},'_filt_','180_','.mat');
wpath=strcat(writepath,'2flt\',flag{1},'\');

im=im2double(fitsread(strcat(rpath,rfile))');
% crop away area of bad DC image
% border pixel calculation
im=im(81:end,201:2250);

% filter gammaspots(median)
im=f_MediFilter1(im,thresh_medi1);

%filter for bright & dark spots
im=f_MediFilter2( im,kernel_medi2,thresh_medi2,'median' );
im(6:end-5,6:end-5); %cout out filter's bad area        


% save
if ind_write==1
    s2_parforSave(strcat(wpath,wfile),im)
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

