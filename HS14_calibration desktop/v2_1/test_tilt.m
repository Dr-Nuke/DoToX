% test fo image tilt
load('xls_file')
flag={'emp';'d2o';'cl3'};
lagging=zeros(3,375);
figure(1);cla;hold on

for i=1:3
    rpath=strcat(writepath,'4cor\',flag{i},'\');
    disp(i)
    for j=1:375
        
        rfile=strcat(flag{i},'_cor_',sprintf('%04d',j),'.mat');
        im=f2_boLoad(strcat(rpath,rfile));
        lower=mean(im(:,196:205),2);
        upper=mean(im(:,1796:1805),2);

        [c,lags]=xcov(lower,upper);
        [~,I] = max(abs(c));
        
        %subpixeL:
        
        subshift=0.5*(log(c(I-1))-log(c(I+1)))/(log(c(I-1))+log(c(I+1))-2*log(c(I)));
        
        lagging(i,j)=lags(I)+subshift;
        
        
        
        
        
        
        
    end
    plot(lagging(i,:))
end
legend('emp','d2o','emp')
title('tilt vs angle step')


; %1

; %2

; %3

; %4

; %5

; %6

; %7

