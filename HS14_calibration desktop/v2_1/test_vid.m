close all
fig=figure(1);
fig.Position=[-1920 0 1920 1080]
cla
rpath='C:\data\tomo_HS14\processed\4cor\cl3\';

for i=1:375

        rfile=strcat('cl3_cor_',sprintf('%04d',i),'.mat');
        load(strcat(rpath,rfile))
        imbo4(var)
        title(num2str(i))
        pause(0.05)
end
    
    
    
    
    
; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %