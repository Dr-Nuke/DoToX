% we rest if we can find the intersection of a circle and a beam


P=f4_ChannelPhan2;

f42_PlotPhan(P,2);

Lx=[0,10];
Ly=[0,0];
line(Lx,Ly,'Color','green')



geom.d=0;               % detector height in mm
geom.nd=1;              % detector pixels
geom.l=50;              % axis-detector distance in mm
geom.Lx=100;            % source-axis distance in mm
geom.Ly=0;              % vertical source offset
geom.nangle=1;          % # of projection angles
geom.noiselvl=0.3;      % noise lvl

geom.angles=linspace(0,360,geom.nangle+1)'; % list of angles
geom.angles(end)=[];

geom.Aheight=geom.nd*geom.nangle;    % height of Matrix A: rotation angles times detector pixels

detxy=f4_DetCoords(geom.d,geom.nd,geom.l);
srcxy=[-geom.Lx;geom.Ly];

figid=1;

%% pre allocate
b=zeros(1,geom.Aheigh,'single'); % phantom imaging data, must be calculated

 
for i=1:length(geom.angles)
    Prot=f42_RotPhan(P,a);
    for j= 1:geom.nd
        for k=1:length(PRot.n)
            
            % find intersects and the two materials
            if P.CL(i)==1 %if the secment is an arc
                [t,m]=f42_IntersecArcLine(P,i,L1(1),L1(2),L2(1),L2(2));
                % t is the beam length [0,1] of an intersect
                % m is the two materials at the intersect
            
            
            else % if the segment is a line    
                [t,m]=f42_IntersecArcLine(P,i,L1(1),L1(2),L2(1),L2(2));
            end
        end
        % after having found all intersects of the beam, we can sort them
        % in ascending t and then find the attenuation langths of the beam
        
        % sort()
        
        % compute individual section lengths
        
        % multiply the section lengths with the attenuation
        
        % add up to find the total attenuation of this (pixel,angle)
        
    end
    % we now have a full row of simulated imaging data
    % do this multiple times to get the sinogram
end

            
            
            
            