function [phan] = f4_ChannelPhan(n)
% creates the channel phantom. based on the 4-symmetry, 
% i.e. computes the values for the first quadrant and rotates other values
% accordingly
% n = amount of segments

t=linspace(0,1,n+1); % sampling parameter

%% geometric data from inventor CAD
r=[13.604,2,5.14,2,13.604];                             % radii
ai=[7.07,125.34,-174.83,125.34,7.07];ai=ai/360*(2*pi);  % angle increment in rad
a0=[0,7.07,-47.59,-42.41,90-7.07];a0=a0/360*(2*pi);     % statr angles in rad
cxy=[0,11.516,6.7,1.429,0;...                           % center points
    0,1.429,6.7,11.516,0];

l=abs(r.*ai);               %length of the segments
lsum=zeros(1,length(l)+1);  %summed length, =1/4 of circumference
for i=1:length(r) 
    lsum(i+1)=sum(l(1:i));
end

%scaling & cropping to 1st quadrant
t=4*t;
t_int=floor(t); % used for rotation later
t=(t-t_int)*lsum(end);

%% attenuation coefficient distribution
att=[0.2,0.5,1,0.5,0.2];

%preallocation
x=zeros(1,n);
y=zeros(1,n);
phan.xy=[x;y];
phan.c=zeros(1,n-1);
%p=zeros(length(t),length(lsum));

for i=1:length(t) %iterate the t parameter, i.e. the segment start- / endpoints
   % if t(i)~=0 %exclude zero index
        p(i,:)=lsum<=t(i); %check which cornerstons come before the curent t
        [pmax(i),pind(i)]=max(lsum(p(i,:))); % check the biggest cornerston amonge them
        lremain(i)=t(i)-lsum(pind(i)); % check how much further t is from the last cornerstone
        aremain(i)=lremain(i)/r(pind(i))*sign(ai(pind(i))); %remaining angle 

        % calculate t's new x and y coordinates
        x(i)=cxy(1,pind(i))+r(pind(i))*cos(a0(pind(i))+aremain(i));
        y(i)=cxy(2,pind(i))+r(pind(i))*sin(a0(pind(i))+aremain(i));
        phan.c(i)=att(pind(i));
        

    %rotate if not in 1st quadrant
    xy=f4_RotMat(t_int(i)*90)*[x(i);y(i)];
    x(i)=xy(1);
    y(i)=xy(2);
end

phan.xy=[x;y];
phan.c(end)=[];
phan.c=rand(1,length(phan.c))+1;


%compute segment length
for i=1:length(t)-1
    phan.l(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
    
end


%phan.c=ones(1,n-1);


if 0 % want plot?
    figure(98)
    subplot(2,3,[1:6])


    plot3(phan.xy(1,1:end-1),phan.xy(2,1:end-1),phan.c,'-x','DisplayName','phantom')
    hold on


    grid on
    xlabel('x')
    ylabel('y')
    zlabel('local attenuation')
    legend('show')
    title(sprintf('%d object segments',n))

end


end

