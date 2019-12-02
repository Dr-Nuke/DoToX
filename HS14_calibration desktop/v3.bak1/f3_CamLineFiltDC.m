function [out] = f3_CamLineFiltDC(im)
% filters out the camera- specific lines
% attention: hardcoded values!!!

% the typical absolute difference of the mean value for a faultiy line is 4500 

% find the lines, summing up the columns
prof=sum(im(1:340,:),1)/340;

%find outliers by means of 3-pix-median
lines=find(abs(prof-medfilt1(prof,3))>3000);

%find the exact line end
for jj =[1:length(lines)] %iterate the lines
    j=lines(jj);
    i=590; %x-ccoordinate
    k=0;
    while k<3 % iterate pixels.
        % if 3 pixels in a row are good, the line ended
        if im(i,j+1)-2*im(i,j)+im(i,j-1)<6000 %if too small difference, check next
            k=k+1;
            i=i+1;
        else
            k=0;
            i=i+1;
        end
    end
    %correct image
    im(1:i-4,j)=(im(1:i-4,j+1)+im(1:i-4,j-1))/2;
    
end
out=im;



end

