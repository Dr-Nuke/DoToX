% here we test the sub pixel shift

%make a phantom

n1=9
n2=10
a=ones(9,10);
[x,y]=ndgrid(1:n1,1:n2);
a(mod(x+y,2)==0)=0;
figure(20);clf;
imshow(a',[0,1]);set(gca,'YDir','normal');

figure(21);clf;


s=linspace(-10,10,11)/10;
ls=length(s);

atest=zeros(n1,n2,ls,ls); % preallocate
for i =6:6 % iterate x shift
    for j=1:ls%iterate y shift
        k=(i-1)*ls+j;
        k=j;
        subplot(1,ls,k);
        shift=[s(i),s(j)];
        shift=[s(j),0];
        %aa=
        atest(:,:,i,j)=interpn(x,y,a,x-shift(1),y-shift(2),'cubic');
        imshow(squeeze(atest(:,:,i,j))',[0,1]);set(gca,'YDir','normal');
        title(sprintf('%.1f %.1f',shift(1),shift(2)))
        ;
        
    end
end
axis tight
        