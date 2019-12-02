%% m mobile

% p = prepay
% m= mini one
p.mb=0.28; %CHF per MB
p.pm=0.28; %CHF per minute
p.mbf=10; % MB free per month

m.mb=0.1; %CHF per MB
m.pm=0.25; %CHF per minute
m.mbf=600; % MB free per month
m.mf=60 % minutes free
m.g=19.90; % grundgebühr


mins=[0:5:120];
MBs=[0:10:50,100:50:600];

[mi,mb]=ndgrid(mins,MBs);

z=zeros(size(mi));
prep=mi*p.pm+max(z,mb-p.mbf)*p.mb;
mo=max(z,mi-m.mf)*m.pm+max(z,mb-m.mbf)*p.mb+m.g;

% MBs surf plot
figure(1);clf;
surf(mi,mb,mb)
xlabel('mins')
ylabel('MBs')
zlabel('MBs')
title('Megabyte plot')

% Mins surf plot
figure(2);clf;
surf(mi,mb,mi)
xlabel('mins')
ylabel('MBs')
zlabel('MBs')
title('Minutes plot')

%prepay cost plot
figure(3);clf;
surf(mi,mb,prep)
xlabel('mins')
ylabel('MBs')
zlabel('CHF')
title('prepay cost')

%mini one cost plot
figure(4);clf;
surf(mi,mb,mo)
xlabel('mins')
ylabel('MBs')
zlabel('CHF')
title('mini one cost')

%difference plot
figure(5);clf;
surfc(mi,mb,prep-mo)
xlabel('mins')
ylabel('MBs')
zlabel('CHF')
title('cot difference prep-mo')
hold on
surf(mi,mb,z)































