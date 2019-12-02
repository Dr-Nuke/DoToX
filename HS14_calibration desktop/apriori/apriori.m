%% this script is my attempt to use apriori information for reconstructing 
% the liquid film thicknesses

% basically we solve Ax=b via ATAx=ATb, (least squares fit) with
% A a matrix containing intersection angles of a beam
% x the unknown wall thicknesses
% b the imaging data
clear all 
clc
format compact

%turn off singular matrix warnings
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:singularMatrix')

% add paths of older functions
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v2_1')
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop/v3')

%% settings
% the geometric (0,0) point is the rotation axis of the object

%geometric settings

geom.d=40;            % detector height in mm
geom.nd=301;          % detector pixels
geom.l=50;            % axis-detector distance in mm
geom.Lx=100;          % source-axis distance in mm
geom.Ly=0;            % vertical source offset
geom.gamma=150;       % # of projection angles

cosphi_thresh=0.3;    % cos phi theshold to consider a beam

%% create detector & source coordinates
detxy=f4_DetCoords(geom.d,geom.nd,geom.l);
srcxy=[-geom.Lx;geom.Ly];

noiselvl=0.3;
figid=1;

%% creat a phantom
% phantoms must specify end points of each section.
% #sections is one less than point pairs, so is attenuation values
%phan=f4_CircPhan(10,1);
phan=f4_ChannelPhan(120);
%phan=f4_SquarePhan(50);
%phan=f4_SimplePhan;

for ii =1:1 % cases iteration
    sprintf('running case %d',i)

    gammas=linspace(0,360,geom.gamma+1)'; % list of angles
    gammas(end)=[];
    %% pre allocate

    A.height=geom.nd*geom.gamma;    % height of Matrix A: rotation angles times detector pixels
    A.width=length(phan.c);         % width of Matrix A: # of object segments
    A.A=zeros(A.height,A.width,'single');   % preallocation

    b=zeros(1,A.height,'single'); % phantom imaging data, must be calculated
    
for i = 1:geom.gamma %iterate all angles
    phan_i=f4_PhanRotate(phan,gammas(i));
    A.pathlengths=0;
    
    for j=1:geom.nd % iterate detector pixesl/ beams
        jj=(i-1)*geom.nd+j;  % counter [1,...,A.height]
        
        for k=1:length(phan_i.c)    % iterate the object sections
                                    % note that is one less than #xy-pairs
                                    % kk=mod(k,length(phan_i.c))+1;
                                     %disp([i,j,k,jj,kk]);
            [t,cosphi]=f4_intercept(srcxy,...
                                    detxy(:,j),...           
                                    phan_i.xy(:,k),...
                                    phan_i.xy(:,k+1));
                                %disp(t');
            %check if t are in [0,1] and if cosphi is good 
            % 
            if all([t>0;t<1;cosphi>cosphi_thresh])
                
                % add the 1/cos to the matrix element  of the [rotation &
                % beam] and [section]
                %A.A(jj,k)=1/cosphi;
                A.A(jj,k)=1/cosphi;
                
                %simulated projection
                b(jj)=b(jj)+(phan_i.c(k)/cosphi)*(1+noiselvl*(2*rand-1));
                % add simulated projection data
            end
        end
    end
    
end
    %%

x=(A.A'*A.A)\((A.A'*b'));

xtot(ii,:)=x;

end
%%
f4_PlotResult(phan,x,figid,geom,noiselvl)