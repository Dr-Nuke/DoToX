%% THis scrip uses the reconstructinos and reconstructs the Liquif Film thickness

%% s3_LFT

if  any(itr_stp==8) % execute this script?
    i=8; % specifies step index of this file, is two
    t1=toc;
 
    [~,name,~]=fileparts(mfilename('fullpath')); %get m file name
    disp(['deploying ',name,'...']) %command line output    
       
    % loard recons
    block_old=strcat(nam_block,nam_stp{i-1},'_',pfx_cas{3});
    path_block=strcat(dir_work,block_old);

    block=imresize(f3_loadSingleVariableMATFile(path_block),refac);

    % F means Film, varia
    F.dr=5;                   % deviation around the radius in pixel
    F.r_min=round(F.r_rod*res-0.5*F.dr);  %physical radius in pixel, with tolerance
    F.r_min=112;               % empirical check shows
    F.S=0.99;                  %sensitivity threshold
    F.o =40;                  % overlap

    F.n_path=round((90/360)*2*pi*F.r_rod*res);  % number of sample paths per pin 
    F.r_path=round(6.7*res)+5; % path length in pixel, res=22.1
   
    F.n_pins = 4; % number of fuel pins
    F.d_angle=135; % angel range
    F.n_angles=F.d_angle+1; % number of angles
    F.h=imrange(2,2); % stack height
    
    
     %preallocate
    c=zeros(imrange(2,2),4,2);  % center [x,y] coordinates of circles
    rad=zeros(imrange(2,2),4);  % radius values of the circles
    m=zeros(imrange(2,2),4);    % confidencey of this circle fit
    F.a0=zeros(imrange(2,2),4); % first angle for each circle
    %a00=zeros(1,4);             % parfor-slicing helper a0
    F.p_end=zeros(imrange(2,2),4);   
    
    centeringrange=f_PixReCalc(500,refac):f_PixReCalc(2100,refac);
    
    % parfor- workaround
    
    
    disp('finding pin centers...')
    parfor k=1:imrange(2,2) %iterate planes
        f_BoCount(k,20,10,5)
        %slice data
        %block(:,:,k)=imrotate(block(:,:,k),15,'crop');
        im=squeeze(block(:,:,k));   
        
        
        %create black/white-image (helps the algorhythm)
        imbw=im2bw(f_normalize(im),0.6); %empirically chosen threshold of 0.6 for greylevel to bw coonversion
        
        % find centers & radii
        try
            [c(k,:,:),rad(k,:),m(k,:)]=f_hugh_4_2(imbw,F.r_min,F.dr,F.S,F.o,0.6); %have this file in the same folder
            if ~all(rad(k,:))
                fprintf ('%4d zero! detected\n',k)
            end

        catch
            fprintf('%4d had a problem assigning the circle values\n',k)
        end
        
        
    end % k-loop
    F.c=c;
    F.rad = rad;
    F.m = m ;
    clear c rad m
    fprintf('\n')
    
    %%     determine centers, do the fit
    F=f3_FindRodCenters(F,centeringrange);
    
    %% create the angles & image profile end points
    F=f3_FindProfilePaths(F);

    %% get image profiles %improfile_tester.m
    %preallocate
    ImProfiles=zeros(F.h,F.n_pins,F.n_angles,F.r_path+1);

    % parfor doesnt like structs, so i copy tha values needed
    h=F.h;
    n_pins=F.n_pins;
    n_angles=F.n_angles;
    cfit=F.cfit;
    ProfEndPnt=F.ProfEndPnt;
    r_path=F.r_path+1;

    fprintf('getting ImProfiles...')
    
    parfor i=1:h
        f_BoCount(i,20,10,5)
         for j=1:n_pins
            for k=1:n_angles

                ImProfiles(i,j,k,:)=...
                reshape(improfile(block(:,:,i)',... % profile of the image
                [cfit(i,j,1),ProfEndPnt(i,j,k,1)],... % from start point
                [cfit(i,j,2),ProfEndPnt(i,j,k,2)],... % to end point
                r_path),...       % with these many samples
                [1,1,1,r_path]);    % and play with reshaping abit  
            end
        end
    end
    
    F.ImProfiles=ImProfiles; clear ImProfiles;
    fprintf('done.\n')
    %% now find the liquid film thicknesses

    % 1st approach: sum up the integrals
    F.LFT=sum(F.ImProfiles(:,:,:,:),4);
    savefast('LFT HS14 diff','F');





end



%% debug scripts

%% plot center coordinates & fit & radii
figure(18)
clf
col=hsv(4);
for k=1:4
    plot(F.c(:,k,1),'Color',col(k,:));
    
    hold on
    plot(F.cfit(:,k,1),'k--');
    plot(F.c(:,k,2),'Color',col(k,:));
    plot(F.cfit(:,k,2),'k--');
    
end

figure(19)
clf

for k=1:4
    plot(F.rad(:,k,1),'Color',col(k,:));
    hold on
    
    
end



%%
figure(20)
clf

for k=100:200%:imrange(2,2)
    figure(18)
    poin=scatter(k*ones(1,8), [F.c(k,:,1),F.c(k,:,2)],'kx');
    
    figure(20)
    clf
    imshow(squeeze(block(:,:,k)'),[]);set(gca,'YDir','normal') % image
    hold on
    %scatter(F.c(k,:,1),F.c(k,:,2),'gx') % center points
    scatter(F.c(k,:,1),F.c(k,:,2),'rx')
    plot(F.c(k,:,1),F.c(k,:,2),'b')
    plot(F.cfit(k,:,1),F.cfit(k,:,2),'g')
    title(k)
    waitforbuttonpress;
    figure(18)
    delete(poin)
end
    
    
%%     test f3_FindProfilePaths
for i=100:10:400 % plane
figure(28)
clf
imshow(squeeze(block(:,:,i)'),[]);set(gca,'YDir','normal') % image
hold on
scatter(F.c(i,:,1),F.c(i,:,2),'rx')
plot(F.c(i,:,1),F.c(i,:,2),'b')
plot(F.cfit(i,:,1),F.cfit(i,:,2),'g')
for j=1:4
    for k=1:10:F.n_angles
        plot([F.cfit(i,j,1),F.ProfEndPnt(i,j,1,k)],...
             [F.cfit(i,j,2),F.ProfEndPnt(i,j,2,k)]) %plot beams
    end
end
title(i)
waitforbuttonpress;
end


