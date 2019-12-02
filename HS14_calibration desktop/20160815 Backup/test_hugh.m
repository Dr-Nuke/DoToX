close all
i=550;
    fnn=sprintf('%04d',i); %file name number
    name_fits=strcat('CL3_',fnn,'.fits');
    filepath=strcat(writepath,'CL3\',name_fits);

cl33=im2double(fitsread(filepath)');
cl33=ones(759);

emp=f_normalize(im2double(fitsread('Empty_0040.fits')'));

im=cl33;
%im=edge(cl3,'Canny',0.12).*f_KernelGen(749,690,1);
%im=edge(emp2,'Canny',0.12).*f_KernelGen(749,690,1));
%im=medfilt2(im2bw(emp2,0.2),[5 5]);
k=111;
r=k;
dr=5;
S=0.99;


tic



o =40; %overlap

[c,rad,m]=f_hugh_4(cl33,r,dr,0.9,o,0.6);
toc

disp([c,rad,m])


imbo3(im,1)
hold on
plot(c(:,2),c(:,1),'g')
viscircles(fliplr(c), rad,'EdgeColor','b');
c(5,:)=c(1,:);
plot(c(:,2),c(:,1),'g')


;
