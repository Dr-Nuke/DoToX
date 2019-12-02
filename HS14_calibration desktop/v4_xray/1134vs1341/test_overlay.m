clear all
clc
format compact
close all


% load film movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1134h/imc2.mat')

% load empty movie
load('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v4_xray/1341h/imc1.mat')
%%



a=squeeze(imc1(:,45,:));
b=squeeze(imc2(:,45,1:1455));

figure(1);clf;
subplot(1,2,1)
imshow(a',[]);set(gca,'YDir','normal')
subplot(1,2,2)
imshow(b',[]);set(gca,'YDir','normal')

figure(2);clf;
imshowpair(a,b)
colorbar
%%
amax=size(a);
bmax=size(b);


acropx=1:amax(1);
acropy=1:amax(2);
bcropx=1:bmax(1);
bcropy=1:bmax(2);

figure(3);clf,
imshowpair(a(acropx,acropy),b(bcropx,bcropy))
axis([1,50,420,460])
figure(4);clf;
imshowpair(a(acropx,acropy),b(bcropx,bcropy))



%%
% lets call a the reference and be the to be shifted image
xscovline=425;

[c,lags]=xcov(a(xscovline,:),b(xscovline,:),'coeff');
[~,ind]=max(c);
shift=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));

inte=fix(shift);
frac=shift-fix(shift);
inte2=ceil(abs(shift))*sign(shift); % 
b2=b; % the new image

%shift in frames direction
b2=circshift(b,inte,2);
[XI,YI]=ndgrid(1:size(b2,1),1:size(b2,2));
b3=interpn(XI,YI,b2,XI+double(frac),YI,'cubic'); %subgrid shift



figure(5)
ax(1)=subplot(1,3,1)
imshowpair(a',a');set(gca,'YDir','normal')
ax(2)=subplot(1,3,2)
imshowpair(a',b2');set(gca,'YDir','normal')
ax(3)=subplot(1,3,3)
imshowpair(a',b3');set(gca,'YDir','normal')
linkaxes(ax)


figure(7)
ax(1)=subplot(1,3,1)
imshow(a'-a',[0,0.1]);set(gca,'YDir','normal')
ax(2)=subplot(1,3,2)
imshow(a'-b2',[0,0.1]);set(gca,'YDir','normal')
ax(3)=subplot(1,3,3)
imshow(a'-b3',[0,0.1]);set(gca,'YDir','normal')

linkaxes(ax)
%%

figure(6);clf;
ax=gca();
plot(a(xscovline,:))
hold on
%plot(b(xscovline,:))
plot(b2(xscovline,:))


%% try to brute force xcov
for i=5:amax(1)-5
    disp(i)
    
[c,lags]=xcov(a(i,:),b(i,:),'coeff');
[~,ind]=max(c);
shift(i)=lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
        (log(c(ind-1))-2*log(c(ind))+log(c(ind+1)));
end
figure(21);clf;
plot(shift)




