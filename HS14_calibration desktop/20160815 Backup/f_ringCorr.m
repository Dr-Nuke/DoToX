function [out] = f_ringCorr(in,val)
% does the ring artifact correction
% in = input 2d- or 3d array
% val = smooth value for the smooth function 


if length(size(in))==3 %its the image stack...
    
    m=mean(in,3);
    s=zeros(size(m));
    for i = 1:size(in,2)
        s(:,i)=smooth(m(:,i),val,'loess');
    end
    s=s./m;
    out=repmat(s,1,1,size(in,3));
    out=out.*in;
    
    
elseif length(size(in))==2 % ...or only an image
    
    m=mean(in,2);
    s=zeros(size(m));
    for i = 1:size(in,2)
        s=smooth(m,val,'loess');
    end
    s=s./m;
    out=repmat(s,1,size(in,2));
    out=out.*in;
else
    disp('bad input dimension number')
end



end

