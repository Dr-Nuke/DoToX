
subshift=zeros(2450,1);

ssmin=480;
ssmax=2200;

for l=ssmin:ssmax
    %disp(l)
    subshift(l)=f3_subpixcorr(cl3(:,l,1),cl3(:,l,end));
    
    
end


plot(ssmin:ssmax,subshift(ssmin:ssmax))

title('centering vs. pixel height')
ylabel('axis shift')
xlabel('y-pixel')