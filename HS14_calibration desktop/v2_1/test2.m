
rpathC='C:\data\tomo_HS14\processed\5sin\cl3\'
rpathD='C:\data\tomo_HS14\processed\5sin\d2o\'
flag={'emp';'d2o';'cl3'};


for j=2000:2050
    rfileC=strcat(flag{3},'_sin_',sprintf('%04d',j),'.mat');
    load(strcat(rpathC,rfileC));
    cl3=var;
    
    rfileD=strcat(flag{2},'_sin_',sprintf('%04d',j),'.mat');
    load(strcat(rpathD,rfileD));
    d2o=var;
    
    if sum(sum(abs(imag(-log(cl3./d2o))))) ~=0
        disp(j)
    end
    
%     imbo3((-log(cl3./d2o))',1);caxis([0,2])
%     title(sprintf('%i',j));
%     pause(0.001)
end
; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

; %10

