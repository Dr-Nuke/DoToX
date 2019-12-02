%% this script is intended to do the a priori tomography, with film
% 1) create a phantom
% 2) create a sinogram from it
% 3) compare it with a measured sinogram
% 4) iteratively find center location (x,y), p0 (rotation angle), detector
%       distance and angle, and source distance

clear MuS P P2 d t L
tic
t1=tic;
P=f42_ChannelPhanTransform_w_film;
P=f42_RotPhan(P,deg2rad(5)); % prevents most singular x= A\b warnings
%P=f42_PhanMove(P,0,30*0.045);
srcxy=[-1000,0];

det_y =31; % width of detector in mm
ndet=251; % number of detector pixels
detpix=linspace(-det_y/2,det_y/2,ndet);

nrot=376;
angles=linspace(0,2*pi,nrot+1);
angles(end)=[];

% my values for cross sections
% region swich: 1:air 2:AL 3:D2 4: inside (ChCl3vapor)
%in 1/mm, for xray
xs=[0,...           % air
    0.09212,...      % Alu
    0.02075,....     % H20
    0.07708,...     % Cl3 liq.
    0.00024];     % Cl3 vapor

c=[[0.7,1,0.7]; % line colors
    [0.4,0.4,0.4];
    [0,0,1];
    [1,0,0.5]];

MuS=zeros(nrot,ndet);

d=zeros(ndet*nrot,5);
%pre-calculate beams     

for i=1:ndet
    eend=[25,detpix(i)];
    l = f42_CreateLineObj(srcxy,eend);
    L{i}=l;
end
    
for i=1:nrot %iterate angles
    P2=f42_RotPhan(P,angles(i));
    disp(i);
    for j=1:ndet % iterate beams
        kk=(i-1)*ndet+j;
        %disp([i,j])
        [t,xy,n_seg,regio,s] = f42_PhanLineIntersect(P2,L{j});
        d(kk,:)=f42_RegioDistances(t,regio,L{j});
        MuS(i,j)=exp(-xs(regio)*s); % accumulatedattenuation
    end
    
%     for j=1:size(xy2)-1
%         l=f42_CreateLineObj([xy2(j,1),xy2(j,2)], [xy2(j+1,1),xy2(j+1,2)]);
%         f42_PlotLineObj(l,ax,[1],c(regio(j),:));
%     end
end
    % end


s42_compareDataToSim
f42_SegmentHistogram
figure(11);cla;
f42_PlotPhan(P,gca,[1,0,0],[1])
for i=1:ndet;f42_PlotLineObj(L{i},gca,[1,2],[0,1,0]);end

t2=toc;
disp(t2-t1)