function [col] = f4_Colorgen(cmap,cmin,cmax,x)
%returns a color vector according to given min & max values and a color map



t=round(((x-cmin)/(cmax-cmin))*size(cmap,1));
t=max(t,1);

col=cmap(t,:);

end

