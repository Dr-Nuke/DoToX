function [out] = f_MediCustom(im,kernel,mode)
% filters an image according to given kernel
% mean and median mode allowed

[xi,yi]=size(im); %image size
[xk,yk]=size(kernel); %kernel size

if and(mod(xk,2)==0,mod(yk,2)==0) %check kenel size
    out=im;
    disp('bad kernel size, filter aborted')
    return
end

rk=floor([xk,yk]/2);    % "radius" of kernel  (distance to image boundary)
rk2=ceil([xk,yk]/2);    % same as rk, but counts also the middle point

%preallocation
tmp=zeros([size(im)+[2*rk(1),2*rk(2)],sum(kernel(:)~=0)]);

ind=find(kernel); %number of non-zero entries
[x,y]=find(kernel); % coordinates of non-zero entries

for i = 1:length(ind) %iterate kernel pixel
    tmp(x(i):xi+x(i)-1,y(i):yi+y(i)-1,i)=im; % add the shifted image to stack
end

if strcmp('mean',mode) % mean & median
    out=mean(tmp,3);
    
elseif strcmp('median',mode)
    out=median(tmp,3);
else
    disp('bad filter string')
    out=im;
    return
end
out=out(rk2(1):end-(rk2(1)-1),rk2(2):end-(rk2(2)-1));  % crop to original size 

end
    
%% down there the old version
% if strcmp('mean',mode)
%     
%     for i = rk(1)+1:xi-rk(1) %iterate xpix
%         
%       for j =rk(2)+1:yi-rk(2)   %iterate ypix
%           a=im(i-rk(1):i+rk(1),j-rk(2):j+rk(2)).*kernel;
%           out(i,j)=mean(a(kernel~=0));
%       end
%     end
% elseif strcmp('median',mode)
% 
%     for i = rk(1)+1:xi-rk(1) %iterate xpix
%         for j =[rk(2)+1:yi-rk(2)]   %iterate ypix
%          a=im(i-rk(1):i+rk(1),j-rk(2):j+rk(2)).*kernel;
%            out(i,j)=median(a(kernel~=0));
%         end
%         d