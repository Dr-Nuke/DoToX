function f42_SegmentHistogram(d)
% visualizes the histogram distribution

figure(10);clf

nhist=1000; % 

col=[0.5,0.5,0.5;... %alu
    0,0,1;...   %h2o
    1,0,0;...   % cl3 liq
    1,0,0.5];   % cl3 vap

name={'Alu','H2O','CHCl3 liq.','CHCl3 vap.'};
for i=1:4
    subplot(2,2,i)
    if i==3
        bins=linspace(0.38,0.5,nhist+1);
        h(i)=histogram(d(d(:,i+1)>0,i+1),bins,'facecolor',col(i,:));
    else
        
        h(i)=histogram(d(d(:,i+1)>0,i+1),nhist,'facecolor',col(i,:));
    end
    grid on
    if i>2
    xlabel('segment length [mm]')
    end
    
    if  mod(i,2)==1
        ylabel('counts')
    end
    
    legend(name(i))
        
end

currentFigure = gcf;
title(currentFigure.Children(end), 'lengment lengt distributions for a tomogram');

end

