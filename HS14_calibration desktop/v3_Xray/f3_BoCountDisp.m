function [] = f3_BoCountDisp(i,j)
%F3_BOCOUNTDISP prints i to the command line if it is a multiple of j
    if mod(i,j)==0 %command line output
         fprintf('_%d ',i)
         
    end

end

