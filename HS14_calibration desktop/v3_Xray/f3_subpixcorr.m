function [subshift] = f3_subpixcorr(s1,s2)
    %F3_SUBPIXCORR does the sub-pixel correction for the cross correlation
    % typically used with the tilt correction and centering for tomogaphies
    %s1 = series 1, s2 = series 2

    [c,lags]=xcov(s1,s2,length(s1)/2,'coeff');
    
    
    [~,ind]=max(c);
    lag_max=lags(ind);

    %subgrid shift
    subshift=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/(log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
    


end

