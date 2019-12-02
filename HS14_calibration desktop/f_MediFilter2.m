function [ im_out] = f_MediFilter2( im_in,kernel, threshold,mode )
% median filter for intermediate sized bright and dark spots (up to 8 pixel
% wide)

im_back=im_in;
%apply filter
filt=f_MediCustom(im_in,kernel,mode);

% remove unwanted filtering
% 1) map defects
diffbw=abs(im_in-filt)>threshold*im_back;
CC = bwconncomp(diffbw, 4);
L = labelmatrix(CC);

% filter by area size
areamin=6;
S = regionprops(CC, 'Area');
smallcut = ismember(L, find( [S.Area] >= areamin));

%filter by eccntricity
ecc=0.8;
S = regionprops(CC, 'eccentricity');
ecccut = ismember(L,find([S.Eccentricity] < ecc)); 

% combine filters
filter_mask=ecccut.*smallcut;

% correct output
im_in(filter_mask~=0)=filt(filter_mask~=0);

im_out=im_in;

end
%%




%% old version

% function [ im_out] = f_medifilter1( im_in,kernel, threshold,mode )
% % median filter for intermediate sized bright and dark spots (up to 8 pixel
% % wide)
% 
% 
% %apply filter
% filt=f_MediCustom(im_in,kernel,mode);
% 
% % make the difference
% diff=im_in-filt;
% 
% % correct output
% im_in(abs(diff)>threshold*im_in)=filt(abs(diff)>threshold*im_in);
% 
% im_out=im_in;
% 
% end


%1
%2
%3
%4
%5
%6
%7
%8
%9
%10
%11
%12
%13 13
