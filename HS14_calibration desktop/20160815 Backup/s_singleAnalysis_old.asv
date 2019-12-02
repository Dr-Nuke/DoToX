% we now operate with a struct
close all
clear A
clc
%%
A.im=f_normalize(im2double(fitsread('CL3_1940.fits')'));
A.im(1:50,101:300)=0;




% find center
r=112;  % minim radius to look for circular structures
dr=5;   % radius range to look for circular structures
S=0.99; % picky-ness of the algorythm
o=40;   % quadrant-overlap (when the quadrants are not exactly located)

thresh_bw=0.6;  % black-white conversion threshold
A.bw = imcomplement(im2bw(A.im,thresh_bw));

% get the 4 circles & radii
[A.cs.center,A.cs.radii,A.cs.metric]=f_hugh_4(A.bw',r,dr,S,o);

% creat paths / find path points
figure(1); imbo3(A.im,1); hold on
A.paths.da=-90;  % width of angle range in deg
A.paths.n=10;  % number of sample path in #
A.paths.r=150; % path length in pixel

% get the paths
for i=1:4,
    j=mod(i,4)+1;   
    % the starting angle is the direction from one center to the next one
    A.paths.startangle{i}=radtodeg(atan2(A.cs.center(j,2)-A.cs.center(i,2),...
        A.cs.center(j,1)-A.cs.center(i,1)));
    
    % generate the path endpoints
    A.paths.endpoints{i} = f_LinGen(A.cs.center(i,:), A.paths.startangle{i},...
        A.paths.da, A.paths.n,A.paths.r);
    
    % debug-plot
    for k=1:size(A.paths.endpoints{i},1),
        plot([A.cs.center(i,1),A.paths.endpoints{i}(k,1)],...
            [A.cs.center(i,2),A.paths.endpoints{i}(k,2)])
    end
end

% improfile
mult=1; % multiplied with path length will be # of sampling points
for i =1:4 % for 4 centers
    for j=1:A.paths.n %for all endpoints
        [cx,cy,A.profiles.profile(i,j,:)]=improfile(A.im',...
            [A.cs.center(i,1),A.paths.endpoints{i}(j,1)],...
            [A.cs.center(i,2),A.paths.endpoints{i}(j,2)],...
            A.paths.r*mult);
    end
end
dist=sqrt((cx-A.cs.center(4,1)).^2+(cy-A.cs.center(4,2)).^2)/res;

% deduce thickness
r_wall=5.14; % radius of outter wall surface

for i=1:4 % centers
    for j=1:A.paths.n % angles
        
        % find walls end
        [c,A.profiles.idx_wall]=min(abs(dist-r_wall));
        % find the corresponding profile value
        val=A.profiles.profile(i,j,A.profiles.idx_wall);
        %find last profile indes whos values is similar
        [c2,A.profiles.idx_film(i,j)]=find(A.profiles.profile(i,j,:)>...
            A.profiles.profile(i,j,A.profiles.idx_wall),1,'last');
        A.thickness(i,j)=sum(A.profiles.profile(i,j,...
            A.profiles.idx_wall:A.profiles.idx_film(i,j)));
    end
end
    

;

    
    
