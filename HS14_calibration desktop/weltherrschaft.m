figure(1)

xlabel('%')    


for i =1:100
    barh([3+rand()])
xlim([0,100])
% For R2014a and earlier:
ax = gca;
set(ax,'XTick',[])
set(ax,'YTick',[])
title('Weltherrschaft')
xlabel('%')    
    pause(0.05)

end