% function [ output_args ] = f4_ChannelPhan(t)
% % creates the parameter based channel phantom
clear all

n=101;
t=linspace(0,1,n);
%t=0.5;
if any([min(t)<0,max(t)>1])
    error('t is outside of [0,1]')
end


r=[13.604,2,5.14,2,13.604];             % radii
ai=[7.07,125.34,-174.83,125.34,7.07];ai=ai/360*(2*pi);   % angle increment in rad
a0=[0,7.07,-47.59,-42.41,90-7.07];a0=a0/360*(2*pi);    % statr angles in rad
cxy=[0,11.516,6.7,1.429,0;...
    0,1.429,6.7,11.516,0];

l=abs(r.*ai);
lsum=zeros(1,length(l)+1)
for i=1:length(r)
    lsum(i+1)=sum(l(1:i));
end

t=4*t;
t_int=floor(t);
t_old=t*lsum(end);
t=(t-t_int)*lsum(end);


start=[13.604,0];
x=zeros(1,n);
y=zeros(1,n);

for i=1:length(t) %iterate the t values

    
    if t(i)~=0
        p(i,:)=lsum<t(i); %check which cornerstons come before the curent t
        [pmax(i),pind(i)]=max(lsum(p(i,:))); % check the biggest cornerston emonge them
        lremain(i)=t(i)-lsum(pind(i)); % check how much further t is from the last cornerstone
        aremain(i)=lremain(i)/r(pind(i))*sign(ai(pind(i)));

        x(i)=cxy(1,pind(i))+r(pind(i))*cos(a0(pind(i))+aremain(i));
        y(i)=cxy(2,pind(i))+r(pind(i))*sin(a0(pind(i))+aremain(i));
        

    else
        
        
        x(i)=start(1);
        y(i)=start(2);
    end
        xy=f4_RotMat(t_int(i)*90)*[x(i);y(i)];
        x(i)=xy(1);
        y(i)=xy(2);
end

    % in which interval are we?
    % get start point and start angle
    % apply remaining path


disp([x',y'])
figure(1)
hold on
plot(x,y,'-x')
grid on
axis equal
% xlim([0,14])
% ylim([0,14])

%%
for k=1:5
    plot(cxy(1,k),cxy(2,k),'o')
end

