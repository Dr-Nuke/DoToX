function [ im] = f_fits( path, thresh, xpix, ypix)
% reads in a fitsfile and applys filters to it

% read in file 
im=im2double(fitsread(path)');
%disp(strcat(path,' ',num2str(f_borange(im))))
%crop

%% debuggung for zeros and nans

if any(isnan(im(:))),
    disp(strcat('NAN occured in',path ))
end

% if any(im(:)<=0),
%     disp(strcat('less or equal zero occured in',path ))
% end
% 
% if any(im(:)==0),
%     disp(strcat('zero occured in',path ))
% end

%% crop & gammaspot filter
im=im(xpix,ypix);

% filter gammaspots(median)
im=f_MediFilter1(im,thresh);

%line filter for special camera lines
im=f_CamLineFilt(im);

%filter for bright & dark spots
kernel=f_KernelGen(9,9,7.5);
im= f_MediFilter2( im,kernel,0.06,'median' );



end

