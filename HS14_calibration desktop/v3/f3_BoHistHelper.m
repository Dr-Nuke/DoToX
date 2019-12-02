function k = f3_BoHistHelper(Values,thresh)
% finds the index of vector Values at which the cumulative sum exceeds 
% thres*sum(Values)
k=1;
kmax=length(Values);
s0=sum(Values);
s=0;
while and(s<thresh*s0, k<=kmax)
    s=s+Values(k);
    k=k+1;
end


end

