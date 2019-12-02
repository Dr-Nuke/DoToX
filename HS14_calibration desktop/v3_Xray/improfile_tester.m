% tests the improfiles


ImProfiles=zeros(F.h,F.n_pins,F.n_angles,F.r_path+1);

% parfor doesnt like structs, so i copy the values needed
h=F.h;
n_pins=F.n_pins;
n_angles=F.n_angles;
cfit=F.cfit;
ProfEndPnt=F.ProfEndPnt;
r_path=F.r_path;

parfor i = 1:h % iterate images
    for j=1:n_pins % iterate pins
        for k=1:n_angles % iterate angles
            
            ImProfiles(i,j,k,:)=...
                reshape(improfile(block(:,:,i),... % profile of the image
                [cfit(i,j,1),ProfEndPnt(i,j,k,1)],... % from start point
                [cfit(i,j,2),ProfEndPnt(i,j,k,2)],... % to end point
                r_path+1),...       % with these many samples
                [1,1,1,r_path+1]);    % and play with reshaping abit                          
        end
    end
end

%% backup
F.ImProfiles=zeros(F.h,F.n_pins,F.n_angles,F.r_path+1);


a1=F.h
a2=F.n_pins
a3=F.n_angles


parfor i = 1:a1 % iterate images
    f_BoCount(i,1,10,5);
    for j=1:a2 % iterate pins
        for k=a3 % iterate angles
            %F.ImProfiles(i,j,k,:)=...
            a=...
                reshape(improfile(block(:,:,i),... % profile of the image
                [F.cfit(i,j,1),F.ProfEndPnt(i,j,k,1)],... % from start point
                [F.cfit(i,j,2),F.ProfEndPnt(i,j,k,2)],... % to end point
                F.r_path+1),...       % with these many samples
                [1,1,1,F.r_path+1]);    % and play with reshaping abit                          
            F.ImProfiles(i,j,k,:)=a;

        end
    end
end





%%
x=10; % x-dimension
y=11; % y-dimension  
z=12; % z-dimension
L=5;  % sample length
a=zeros(x,y,z,L); % preallocation

parfor i=1:F.h
    for j=1:F.n_angles
        for k=1:z
            disp([i,j,k]) %debug display
            
            a(i,j,k,:)=...
            reshape(rand(L,1),...
            [1,1,1,L]);

        end
    end
end
            
            