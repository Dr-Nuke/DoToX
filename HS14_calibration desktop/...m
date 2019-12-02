function [ out] = f_hugh(im,kernel)
%manual attempt of hugh transform


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
out=zeros([size(im)+[2*rk(1),2*rk(2)]]);

ind=find(kernel); %number of non-zero entries
[x,y]=find(kernel); % coordinates of non-zero entries

for i = 1:length(ind) %iterate kernel pixel
    
    out(x(i):xi+x(i)-1,y(i):yi+y(i)-1)=out(x(i):xi+x(i)-1,y(i):yi+y(i)-1)+im; % add the shifted image to stack
end

out=out(rk2(1):end-(rk2(1)-1),rk2(2):end-(rk2(2)-1))/length(ind);  % crop to original size 




end

