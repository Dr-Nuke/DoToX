xx=[0 5000 10000 15000 20000 30000 50000 100000 150000 200000 250000];;
yy=[25 35 50 70 95 130 180 270 350 400 500];
size(x);
size(y);

x=logspace(log10(1),log10(xx(end)),1000);

y=interp1(xx,yy,x,'previous');
loglog(x,y./x)

grid on



;
