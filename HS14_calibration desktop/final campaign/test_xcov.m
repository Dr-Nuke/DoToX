pixrow=200;
i1=1;
i2=4;

figure(1);clf;
imshowpair(T.sino.Raw{i1},T.sino.Raw{i2})

figure(2);clf
ax=gca();
hold on

l1=T.sino.Raw{i1}(pixrow,:);
l2=T.sino.Raw{i2}(pixrow,:);

fprintf('sizes are %d and %d\n',size(l1,2),size(l2,2))

plot(l1,'Displayname','Reference')
plot(l2,'DisplayName','to be shifted by...')
grid on
legend
maxlag=50;

[c,lags]=xcov(l1,l2,maxlag,'coeff');
if any(isnan(c)) %debug
    fprintf('%4d contains NaN\n',i)
    
else % normal case
    fprintf('%4d does not contain NaN\n',i)
    [~,ind]=max(c);
    shift=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
end

title(sprintf('shift = %f',shift))

%%
[d,T]=f.ReadXrayTomo(T);