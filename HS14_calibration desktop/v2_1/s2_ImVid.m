
flag={'emp';'d2o';'cl3'};
figure(1)
cla
i=3
for j =1:375
    rfile=strcat(flag{i},'_add_',sprintf('%04d',j),'.mat');
    rpath=strcat(writepath,'3add\',flag{i},'\');
    imbo3(s2_boLoad(strcat(rpath,rfile)),1);
    pause(0.05)
end
;
;
;