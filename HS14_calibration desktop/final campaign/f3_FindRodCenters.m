function F = f3_FindRodCenters(F,range)
% fits the circle finder's data to a linear function
fprintf('fitting pin centers... ') 

for i=1:F.n_pins %4 pins
    for j=1:2 %2 center coordinates
        F.p{i,j}=polyfit(range,squeeze(F.c(range,i,j))',1); %linear fit
        F.cfit(:,i,j)=polyval(F.p{i,j},1:size(F.c,1));
    end
end

fprintf('done.\n')
end

