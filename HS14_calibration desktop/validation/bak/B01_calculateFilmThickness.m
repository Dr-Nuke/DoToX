
function [filmThicknessData,integralValues] = B01_calculateFilmThickness(I)

%close all
%clear
%clc

% load ./tomoData/img_empty.mat
% empty = tomoNoisy;
% load ./tomoData/img_film.mat
% film = tomoNoisy;
% I = (film-empty);
%imshow(I,[]);

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
for theta = 0:2:90
    x2 = x+cosd(theta)*L;
    y2 = y+sind(theta)*L;
    %line([x,x2],[y,y2])
end
%figure;hold on
counter = 0;
for theta = 0:2:90
    counter = counter+1;
    x2 = x+cosd(theta)*L;
    y2 = y+sind(theta)*L;
    [a]=improfile(I,[x,x2],[y,y2],L*4);
    profileData(:,counter) = a;
end

integralValues = sum(profileData,1);

profileData = mean(profileData,2);
profileData = profileData(180+1:end);
filmThicknessData = profileData;
%plot(profileData)