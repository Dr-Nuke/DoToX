
figure(1);clf;
x=linspace(0,4000000,1001);
tax=VSS(x);
tax3a=VSS(max(0,x-600000))
taxsave=tax-tax3a;
plot(x,taxsave./x*100)
grid on
title('annual saving of Vermögenssteuer with 600k in pillar 3a')
xlabel('net worth [CHF]')
ylabel('tax savings [%]')


function [tax] = VSS(x)

for i=1:length(x)

v = [0, 77000, 308000, 694000, 1310000, 2235000, 3158000];
s = [0, 0, 115, 501, 1425, 3274, 5580];
add = [0, 0.5, 1, 1.5, 2, 2.5, 3];


t=v-x(i);
t=t(t<=0);
ind(i)=length(t);
V(i)=v(ind(i));
S(i)=s(ind(i));
Add(i)=add(ind(i));
tax(i)=s(ind(i))+(x(i)-v(ind(i)))/1000*Add(i);
end


end

