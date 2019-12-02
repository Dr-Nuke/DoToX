function [centers,radii,metric] = f_hugh_4(im,r,dr,S,o,thresh_bw)
%find the centers of the cooling channels

im=f_normalize(im);
[x,y]=size(im);

% create indices for sub-images that contain a full circle each
x1=[1,ceil(x/2)-o,ceil(x/2)-o,1];
x2=[ceil(x/2)+o,x,x,ceil(x/2)+o];
y1=[1,1,ceil(y/2)-o,ceil(y/2)-o];
y2=[ceil(y/2)+o,ceil(y/2)+o,y,y];

%preallocate
centers=zeros(4,2);
radii=zeros(4,1);
metric=zeros(4,1);




for i=1:4
    S2=S;
    while true % iterate the sensitivity until it is sufficient
    
    im2=im2bw(im(x1(i):x2(i),y1(i):y2(i)),thresh_bw);
    [c,rad,m] = imfindcircles(im2,[r,r+dr],'Sensitivity',S2,'ObjectPolarity','dark');
    
    
        if length(rad)~=0
            centers(i,1)=c(1,2)+x1(i)-1;
            centers(i,2)=c(1,1)+y1(i)-1;
            radii(i)=rad(1);
            metric(i)=m(1);
            %fprintf (' %.4f ', S2)
            break
        elseif S2 > 0.9999
            centers(i,1)=0;
            centers(i,2)=0;
            radii(i)=0;
            metric(i)=0;
            fprintf ('nor circles found with S= %.3f\n', S2)
            break
           
        else
             S2=S2+0.5*(1-S2);
             %disp(sprintf('sensitivity set to %f',S2))
             %fprintf (' %.3f', S2)
        end
   
    end


end
%fprintf ('\n')
end