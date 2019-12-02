function [out] = f_centCrop(sin,lag,xmin,xmax,xdo_1,xdo_2)
% centers the sinogram and crops from both sides
% sin = incoming sinogram
% lag = the lag of xcov(...)

if lag >0
    crp=min((xdo_1-xmin-lag),(xmax-xdo_2));
    sin([1:lag+crp,end-crp:end],:)=[];
    if lag>crp %pad zeros
        n_z=lag-crp; %number of added zeros
        sin=vertcat(sin,zeros(n_z,size(sin,2)));
           
    end

elseif lag<0
    lag=-lag;
    crp=min((xdo_1-xmin),(xmax-xdo_2-lag));
    
    sin([1:crp,(end-crp-lag):end],:)=[];
    if lag>crp %pad zeros
        n_z=lag-crp; %number of added zeros
        sin=vertcat(zeros(n_z,size(sin,2)),sin);
    end


else
    crp=min((xdo_1-xmin),(xmax-xdo_2));
    sin([1:crp,(end-crp):end],:)=[];
end
out=sin();

end

