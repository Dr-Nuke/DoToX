% Program to make a video out of the combined data of a file from certain
% flow conditions. The video can be flat as seen from the top and the pixel
% correspond to single positions, or, the video is 3D flat or 3D polar.
% Select what video type you want with the number below.

%% Tabula Rasa
clear all
close all
clc
tic     % Measure the time to see if everything is working well

%% Initial Values
format short;

% Choose the length of the video in frames
length_video=1000;
frame=10000;    % Frequency of the measurement (for correct time indication)
% Choose video type: 1 = 2D, 2 = 3D flat, 3 = 3D polar, 4 = all three
video_type=3;

dist=2e-3;      % Distance between two electrodes [m]
R=12;           % Radius of the pipe

% To plot the data  correctly  the position of the sensor is needed as well
% as the Degree per Pixel = 360*(ø*Pi/(# Transmitters * distance between
% them))
deg_per_pixel = 360 * (dist * 1000) / (24 * pi);
% Lowest point estimation from mean film thickness ± 0.5 for UpWard /
% DownWard flow
lowest_point_UW = 17.5;
lowest_point_DW = 18.5;
circumference = (0:35) * deg_per_pixel + ...
    (180-(lowest_point_UW-1)*deg_per_pixel);

% Data must be in a folder and be named in the form 
% deg-XX-ul-YYYY-ug-ZZZ_run?.mat
% XX    = angle
% YYYY  = superfacial liquid velocity (underscore is decimal point)
% ZZZ   = superfacial gas velocity (underscore is decimal point)

path='C:\Users\Lukas\Desktop\Multilayer_LFS\AN-UW-1000Hz\Data\';
dirName = [path];                       % folder path to generate file list
files = dir(fullfile(dirName,'*.mat')); % list all '.dat' files  
files = {files.name}';                  % get file names

% Path for the video to be stored
results_path='C:\Users\Lukas\Desktop\Multilayer_LFS\AN-UW-1000Hz\';


for file_nmbr =1:2:numel(files)  % Change file in that folder. 

% Load the file of which the video should be made. To make it obsolete to
% know the variable's name of the loaded file (e.g. comb_run1) , the
% loading is a little more complicated than in the other files.

filename=char(files(file_nmbr));    % Get name of current file
load_file=load([path, filename]);   % Load that file
disp([path, filename])              % Display the paht of the loaded file
% Get name of the variable (e.g. comb_run1)
data_name=char(fieldnames(load_file));  
data=load_file.(data_name);         % Access loaded data
clear load_file                     % Clear loaded struct

% Titles for each file-name format (depending on the air flow)
% title_name=['u_l = 0.', filename(13:14),...
%     ' m/s, u_g = ', filename(19), '.', filename(21), ' m/s'];
title_name=['u_l = 0.', filename(13:14), ' m/s, u_g = ', ...
    filename(19:20), ' m/s, Inclination = ', filename(5:6), '°'];
 %%
% 3D Polar 
if video_type == 3 || video_type == 4
%Idea: Assign every measured point (ny*nx*timefr) catesian coordinates
%(x,y,z) (-> 3*(ny*nx*timefr)), considering the cylindrical geometry of the
%measurement and the water level at each point. Plot the points taking the
%color from the initial measuring matrix surf(x(t=i),y(t=i),z(t=i),data);
%%
% Initialize Matrices for the plot
x_pol=single(zeros(size(data,1),size(data,2),length_video));       
y_pol=single(zeros(size(data,1),size(data,2),length_video));       
z_pol=single(zeros(size(data,1),size(data,2),length_video));     

% Polar coordinates in axial direction 
axial_values=linspace(0,size(data,1)*dist*1000,size(data,1)); 
% Polar coordinates in tangential direction
phi=(pi-circumference/180*pi)/2;    

r=-data+R;                    % polar coordinates in radial direction

for k=1:length_video                % go through time
    for i=1:size(data,2)            % go through angles/phi          
        x_pol(:,i,k)=axial_values;  % axial coordinates are trivial
        phi_curr=phi(i);            % set current angle
        for j=1:size(data,1)        % go through axial direction  
            y_pol(j,i,k)=sin(phi_curr)*r(j,i,k);    % y=sin(phi)*r
            z_pol(j,i,k)=-cos(phi_curr)*r(j,i,k);   % z=cos(phi)*r
            
        end
    end
end
%%
fig=figure();
colormap(hsv)
for i=1:length_video  
    %%
    i=90;
    surf(x_pol(:,:,i),y_pol(:,:,i),z_pol(:,:,i),data(:,:,i));
    %view(-85+10*sin(i/1000*2*pi),30+10*cos(i/1000*2*pi));
    view(85,50);
    xlabel('[mm]');
    ylabel('[mm]');
    zlabel('[mm]');
    caxis([0 1])
    h = colorbar;
    ylabel(h, 'Film Thickness [-]')
    title([title_name, ', t = ', num2str(i/frame,'%0.3f'), 's'])  
%%
    M(i) = getframe(fig);
    
end;
video=VideoWriter([results_path, filename, '3D_polar.mp4'],'MPEG-4');
open(video)
writeVideo(video,M);
close(video)
close(fig)
toc
end

%% 3D flat
if video_type == 2 || video_type == 4
y=linspace(0,size(data,1)*dist*1000,size(data,1));
x=circumference;
[X, Y] = meshgrid(x,y);
fig=figure();
colormap(hsv)
for i=1:length_video              
    surf(X,Y,data(:,:,i));
    view(200,40);
    pbaspect([1 1 0.1])
    xlabel('Circumferential Position [deg]');
    ylabel('[mm]');
    zlabel('[mm]');
    zlim([0 3.7])
    caxis([0 2])
    h = colorbar;
    ylabel(h, 'Film Thickness [-]')
    title([title_name, ', t = ', num2str(i/frame,'%0.3f'), 's'])  

    M(i) = getframe(fig);
    
end;
video=VideoWriter([results_path, filename, '3D_flat.mp4'],'MPEG-4');
open(video)
writeVideo(video,M);
close(video)
close(fig)
toc
end

%% 2D Plot
if video_type == 1 || video_type == 4
fig=figure();
colormap(hsv)
for i=1:length_video
    image(data(:,:,i),'CDataMapping','scaled')
    caxis([0 2])
    xlabel('# Transmitter');
    ylabel('# Receiver');
    h=colorbar;
    ylabel(h, 'Film Thickness [mm]')
    title([title_name, ', t = ', num2str(i/frame,'%0.3f'), 's'])  
    M(i)=getframe(fig);
end

video=VideoWriter([results_path, filename, '2D.mp4'],'MPEG-4');
open(video)
writeVideo(video,M);
close(video)
close(fig)
toc
end

end % End of ForLoop for all files in the folder  

disp('Program terminated correctly');
% Program finish. I hope the videos are satisfactory...

%               ()
%               <M
%    o          <M
%   /| ......  /:M\------------------------------------------------,,,,,,
% (O[]TSCHUKI[]I:K+}=====<{R}>================================------------>
%   \| ^^^^^^  \:W/------------------------------------------------''''''
%    o          <W
%               <W
%               ()