dcfile=strcat('dc','_add_','.mat');
dcpath=strcat(writepath,'3add\','dc_','\');
doseXrange=[1:250,950:1120];
dc=f2_boLoad(strcat(dcpath,dcfile));

for i =1:10

dc2=f_MediFilter1(dc,i/100);

%
imbo3((dc2-dc)~=0,i);
end
%%


for i =1:10

dc3=f_MediFilter2( dc,kernel_medi2,thresh_medi2,'median' );
imbo3((dc2-dc)~=0,i);
%

end


% => DC images are not contaminated
;
;
;
;

;
;
