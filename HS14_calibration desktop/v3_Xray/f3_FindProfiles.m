function F=f3_FindProfiles(F)
% finds the image profiles 

% preallocate

F.ImProfiles=zeros(F.h,F.n_pins,F.n_angles,F.r_path)

for 1 = 1:F.h % iterate image planes
    for j=1:F.n_pins % iterate pins
        for k=1:F.n_angles % iterate angles
            F.ImProfile(i,j,k,:)=improfile()
            
            
        end
    end
end


end

