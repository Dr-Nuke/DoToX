
% film file
pathemp= 'C:\Users\cbolesch\Desktop\HS14_calibration desktop\v4_xray\1341h\reconEmp.mat'
pathfilm='C:\Users\cbolesch\Desktop\HS14_calibration desktop\v4_xray\1134h\reconV1.mat'

emp=load(pathemp);
film=load(pathfilm);
%%
diff=film.recon-emp.recon;

%%
save('diff','diff')