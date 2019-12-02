
function [filmThicknessData,integralValues] = B01_calculateFilmThickness(I)


% center of circle found visually by trial and error
%hold on
%r = 50;
x = 85;
y = 85;
%th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% plot circles and lines to visualize stuff
%h = plot(xunit, yunit);
L = 65;

for theta = 0.1:0.1:90
    x2 = x+cosd(theta)*L;
    y2 = y+sind(theta)*L;
end
%figure;hold on
counter = 0;
for theta = 0.1:0.1:90
    counter = counter+1;
    x2 = x+cosd(theta)*L;
    y2 = y+sind(theta)*L;
    [a]=improfile(I,[x,x2],[y,y2],L*4);
    profileData(:,counter) = a/4;
end

integralValues = sum(profileData,1);
integralValues = imresize(integralValues,[1 30],'box');

profileData = mean(profileData,2);
profileData = profileData(180+1:end);
filmThicknessData = profileData;
%plot(profileData)

% close all
% 
% subplot(1,2,1)
% plot(integralValues)
% ylim([0 inf])
% 
% subplot(1,2,2)
% plot(profileData)
% ylim([0 inf])

