
%%

close all
clear
clc

load ./tomoData/img_empty.mat
empty = tomoNoisy;
load ./tomoData/img_film.mat
film = tomoNoisy;
I = (film-empty);

% gives film thickness over top left quarter rod in steps of 0.05 mm
[filmThicknessData] = B01_calculateFilmThickness(I);

L = max(size(filmThicknessData));
plot([1:L]*0.25,filmThicknessData)
grid on