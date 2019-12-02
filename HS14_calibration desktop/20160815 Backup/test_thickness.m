

clc
cl33=f_normalize(im2double(fitsread('V:\1phd\codes\reconstruction\HS14_calibration\CL3_1940.fits')'));
%V:\1phd\codes\reconstruction\HS14_calibration

cl33(1:50,101:300)=0;
clbw=im2bw(cl33,0.6);
clbw=imcomplement(clbw);
k=round(5.14*res)-3;
r=k;
dr=6;
S=0.99;
o =40; %overlap


%% find centers
[c,rad,m]=f_hugh_4(clbw,r,dr,S,o);
c2=fliplr(c);
disp([c2,rad,m])
imbo3(cl33 ,1)
hold on
plot(c2(:,1),c2(:,2),'g')
viscircles(c2, rad,'EdgeColor','b');

plot(c2(:,1),c2(:,2),'g')
disp(length(rad));
 
%% find angles

hold on
da=90;
n=100;
r2=150;

for i=1:4,
    j=mod(i,4)+1;
    disp([i,j]);
    a_start=radtodeg(atan2(c2(j,2)-c2(i,2),c2(j,1)-c2(i,1)));
    [px,py] = f_LinGen(c2(i,:), a_start,da,n,r2);
    for k=1:n,
        plot([c2(i,1),px(k)],[c2(i,2),py(k)])
    end
end
axis([0 749 0 749]); axis equal



%% do the line profiles
;
% imbo3(cl33,2)
% hold on
% [cx,cy,v,xi,yi] = improfile(cl33',[c2(4,1),points(1,1)'],[c2(4,2),points(1,2)'],1000);

;
;
