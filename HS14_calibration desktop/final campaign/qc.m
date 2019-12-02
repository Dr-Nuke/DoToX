classdef qc %quality check  functions
    

%% function collection file for the final campaign Dec. 201+
methods(Static) %evil function scope hack
    

    %% recon optimization
function [kern,mask,gm]=kern(rec,threshFac,spot)
    % returns the kern of the mask
    
    mask=zeros(size(rec));
    % create the mask by gradient
    [gm,~]=imgradient(rec);    
    
    %binarize
    bm=gm>threshFac;

    % remove inner boundary
    %bm=~imfill(~imfill(bm~=0,spot),spot);
    
   
    [L,n] = bwlabel(bm); %assign a number to each area
    if n<6
        error('less han 6 areas detected')
    end
    for i=1:n
        s(i)=sum(L(:)==i);
    end
    [so,I]=sort(s,'descend');
    
    ind=[1 2 3 3 3 3]; %
    
    for i=1:6
        mask(L==I(i))=ind(i);
    end
        
    [~,n2]=bwlabel(mask);
    if n2~=6
        error('bad mask')
    end
    kern=sum(mask(:)~=0);
    
end    

function [v,g]=VarianceQuality(rec,mask)
    % variance based quality criterion
    [g,~]=imgradient(rec);
    
    s1=sum(mask(:)==2);
    s2=sum(mask(:)==3);
        
    gs1=g(mask==2); % al-air edge
    gs2=g(mask==3); % al-d2o edge
    
    v=(s1*var(gs1)+s2*var(gs2))/(s1+s2);


end

function [v,g]=VarianceQualityOld(rec,mask)
    % variance based quality criterion
    [g,~]=imgradient(rec);
    
    gs=g;
    gs(mask==0)=[];
    v=sum(gs.^2);
end
    
end
end



