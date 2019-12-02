
%F4_PLOTERRORS Summary of this function goes here
%   Detailed explanation goes here
figure(54)
clf
subplot(2,1,1)
for i=1:10
    
dx=(xtot(i,:)-phan.c)./phan.c;    




plot(dx,'DisplayName','|deviation|',...
    'Color',f4_Colorgen(jet,1,10,i),...
    'DisplayName',strcat(num2str(10*i),' pixels'))
hold on
grid on

xlabel('object segment#')
ylabel('relative deviation ')
xlim([1,max(length(x),2)])
title('The influence of projection angles numbers on the error')
ylim([-0.2,0.2])

end

legend('show')

subplot(2,1,2)
plot(10:10:100,std(xtot'-repmat(phan.c,10,1)'))
xlabel('#projection angles')
ylabel({'standard deviation',' of the segement-wise error'})