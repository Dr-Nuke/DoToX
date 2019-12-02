function F=f3_FindProfilePaths(F)
% creates the profile end points 
fprintf('creating improfile coordinates... ')

a=1:F.n_pins; % pin indicator
b=circshift(a,-1); %shifted indicator for pin_n minus pin(n-1)


%starting angles
F.a_start=rad2deg(atan2(F.cfit(:,b,2)-F.cfit(:,a,2),...  %y-difference
                        F.cfit(:,b,1)-F.cfit(:,a,1)))... %x-difference
                        -(F.d_angle-90)/2; % add the more than 90deg section

                    
                    
% the angle map
a=linspace(0,F.d_angle,F.n_angles);% generate incremental angles first

F.ProfAngles= degtorad(repmat(F.a_start,[1,1,F.n_angles])... % start angles
    +repmat(reshape(a,[1,1,F.n_angles]),F.h,4,1)); % plus incremental angles

% x values of end points
F.ProfEndPnt(:,:,:,1) = repmat(F.cfit(:,:,1),[1,1,F.n_angles])... % start at center x
                        +F.r_path*cos(F.ProfAngles);     % and move r_path into Profangles direction
F.ProfEndPnt(:,:,:,2) = repmat(F.cfit(:,:,2),[1,1,F.n_angles])... % start at center x
                        +F.r_path*sin(F.ProfAngles);     % and move r_path
                    
fprintf('done.\n')                    
end
                    
%% old, for loop based approach
%                     
% for i=1:size(F.c,1) % iterate plane in z-direction
%     
%     for j=1:F.n_pins
%         jj=mod(j,F.n_pins)+1; % counter for the next cyclical pin
%         %disp([i,j,jj]); %debug
%         
%         F.a_start(i,j)=rad2deg(atan2(F.cfit(i,jj,2)-F.cfit(i,j,2),...
%                               F.cfit(i,jj,1)-F.cfit(i,j,1)))...
%                               -(F.d_angle-90)/2; %calculate start angle
%     
%         [F.ProfEndPnt(i,j,1,:),F.ProfEndPnt(i,j,2,:)] =...
%             f_LinGen(F.cfit(i,j,:), F.a_start(i,j),...
%             F.d_angle,F.n_angles,F.r_path);
%     end
% end
% 
% end

%%


