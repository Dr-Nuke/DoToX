%% this script is intended to do the a priori tomography
% 1) create a phantom
% 2) create a sinogram from it
% 3) compare it with a measured sinogram
% 4) iteratively find center location (x,y), p0 (rotation angle), detector
%       distance and angle, and source distance

tic
t1=tic;
P=f42_ChannelPhanTransform;
%P=f42_RotPhan(P,deg2rad(5));
%P=f42_PhanMove(P,0,30*0.045);
srcxy=[-1000,0];

ndet=301;
detpix=linspace(-15,15,ndet);

nrot=50;
angles=linspace(0,2*pi,nrot+1);
angles(end)=[];

% my values for cross sections
% region swich: 1:air 2:AL 3:D2 4: inside (ChCl3vapor)
%in 1/mm
xs=[0,...           % air
    10.143,...      % Alu
    45.039,....     % D20
    0.485]/1000;     % Cl3 vapor

c=[[0.7,1,0.7]; % line colors
    [0.4,0.4,0.4];
    [0,0,1];
    [1,0,0.5]];

MuS=zeros(nrot,ndet);

d=zeros(ndet*nrot,4);
%pre-calculate beams     
clear L
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

t2=toc;
disp(t2-t1)