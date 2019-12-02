
close all
clear
clc

counter = 0;
% film case is film thickness * 0.1 mm, available points are 0.1:0.1:1.5 mm, or cases 1:15
for filmCase = 1:5
% which Cu thickness to use for X-ray filtering...
%         case: 1   2   3   4   5   6   7   8
% Cu thickness: 0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0
for filterCase = 1:8

counter = counter+1; 
load(['./tomoData/tomo_filmCase_' num2str(filmCase) '_filterCase_' num2str(filterCase) '_polychromatic_withFilm.mat'],'tomoNoisy');
tomoWithFilm = tomoNoisy;
load(['./tomoData/tomo_filmCase_' num2str(filmCase) '_filterCase_' num2str(filterCase) '_polychromatic_empty.mat'],   'tomoNoisy');
tomoEmpty = tomoNoisy;
tomoDiff = tomoWithFilm - tomoEmpty ;
save(['./tomoData/tomo_filmCase_' num2str(filmCase) '_filterCase_' num2str(filterCase) '_polychromatic_diff.mat'],'tomoDiff')

x{filmCase,filterCase} = tomoDiff;
[filmDataAll{filmCase,filterCase},integralDataAllRaw{filmCase,filterCase}] = B01_calculateFilmThickness(tomoDiff);
%subplot(5,8,counter)
%imshow(tomoDiff,[])

end
end


imgDiffAll = ...
   [x{1,1} x{1,2} x{1,3} x{1,4} x{1,5} x{1,6} x{1,7} x{1,8};...
    x{2,1} x{2,2} x{2,3} x{2,4} x{2,5} x{2,6} x{2,7} x{2,8};...
    x{3,1} x{3,2} x{3,3} x{3,4} x{3,5} x{3,6} x{3,7} x{3,8};...
    x{4,1} x{4,2} x{4,3} x{4,4} x{4,5} x{4,6} x{4,7} x{4,8};...
    x{5,1} x{5,2} x{5,3} x{5,4} x{5,5} x{5,6} x{5,7} x{5,8}];

close all
imshow(imgDiffAll,[])
xlabel('filter case 1-8')
ylabel('film thickness case 0.1-0.5 mm')

%%
figure

CuThicknessList = [0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0];

for j=1:8
subplot(2,4,j)
hold on
grid on
xlabel('radial position [mm]')
ylabel('average attn [1/pixel]')
title([num2str(CuThicknessList(j)) ' mm Cu thickness, thickness [mm] in legend'])
for i=1:5
    plot([1:80]*0.025,filmDataAll{i,j})
    integralDataAll(i,j) = sum(filmDataAll{i,j});
end
legend(num2str([0.1:0.1:0.5]'))
ylim([0 0.025])
end

%%

figure;
cvec = 'rgkbkmkk';
hold on
grid on

%         case: 1   2   3   4   5   6   7   8
% Cu thickness: 0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0
title('dashed = 1:1 , lines = extrapolated from first point , markers = simulation')
for j=[1 2 4 6 8]
m = integralDataAll(1,j)/0.1;
plot([0 0.5],[0 0.5*m],cvec(j),'LineWidth',2);
end
for j=[1 2 4 6 8]
scatter([0.1:0.1:0.5]',integralDataAll(:,j),cvec(j),'LineWidth',2)
end
plot([0 0.5],[0 0.5],'--k')
xlim([0 inf])
ylim([0 inf])
xlabel('phantom film thickness [mm]')
ylabel('calculated film thickness from tomogram [mm]')
legend(num2str((CuThicknessList([1 2 4 6 8]))'))

%%

figure;
for i=1:5
    subplot(1,5,i)
    hold on
    grid on
    xlabel('radial angle')
    ylabel('film attenuation integral [a.u.]')
    title(['film = ' num2str(i/10) ' mm, legend = Cu thickness'])
    for j=1:8
        plot([3:3:90],integralDataAllRaw{i,j})
    end
    legend(num2str([0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0]'))
end




