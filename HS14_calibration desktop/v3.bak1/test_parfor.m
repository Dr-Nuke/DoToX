
x=100;
y=1200;
z=30;
data=zeros(x,y,z);
dir_r='/media/data/cbolesch/Tomo_HS14/processed/1raw/emp/';


parfor i=1:z;
    k=sprintf('%04d',i);
    file_r=strcat('emp_',k,'_.fits');
    filename=strcat(dir_r,file_r);
    disp(filename)
    im=im2single(fitsread(filename));
    data(:,:,i)=im(1:x,1:y);
end