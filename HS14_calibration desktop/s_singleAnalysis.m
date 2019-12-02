% we now operate with a struct
close all
clear G
clear LFT

% struct convention
    % G global, total struct with properties equall to all
    % p(i) plane i or reconsruction

    % find center parameters
    G.rmin=111;  % minim radius to look for circular structures (try 5.14*res)
    G.dr=5;   % radius range to look for circular structures
    G.S=0.99; % picky-ness of the algorythm
    G.overlap=20;   % quadrant-overlap (when the quadrants are not exactly located)
    G.thresh_bw=0.6;  % black-white conversion threshold
    
    % creat paths / find path points parameters
    G.da=-90;  % width of angle range in deg
    G.n_path=round((90/360)*2*pi*5.14*res);  % number of sample paths per pin 
    G.r_path=round(6.7*res); % path length in pixel, res=22.1

    G.n_pins=4;
    ind_debug=1;

%for i=1700:1740 % planes
iii=1;
 for i=[550:1950]
    fprintf (' %d ',i)

    fnn=sprintf('%04d',i); %file name number
    name_fits=strcat('CL3_',fnn,'.fits');
    filepath=strcat(writepath,'CL3\',name_fits);
    

     ii=i;
     i=1;


    G.p(i).im=(im2double(fitsread(filepath)'));
    G.p(i).im(1:50,101:300)=0;

    %figure(1); imbo3(G.p(i).im,1); hold on
 
    % get the 4 circles & radii, in 4 quadrants one each
    [c,r,m]=f_hugh_4(G.p(i).im',G.rmin,G.dr,G.S,G.overlap,G.thresh_bw);
    if ~all(r)
        fprintf (' zero! ')
    end
        
    for j=1:4
        G.p(i).pin(j).cntx=c(j,1);G.p(i).pin(j).cnty=c(j,2); %center x & y
        G.p(i).pin(j).rad=r(j);
        G.p(i).pin(j).met=m(j);
    end



    % get the paths
    for j=1:G.n_pins
        
        jj=mod(j,4)+1;   
        % the starting angle a0 is the direction from one center to the next one
        % atan2( y_c2- y_c1, x_c2 - x_c1)
        G.p(i).pin(j).a0=radtodeg(...
            atan2(G.p(i).pin(jj).cnty-G.p(i).pin(j).cnty,...
                  G.p(i).pin(jj).cntx-G.p(i).pin(j).cntx));
         
        % generate the path endpoints
        [G.p(i).pin(j).endpts.x,G.p(i).pin(j).endpts.y] = f_LinGen(...
            [G.p(i).pin(j).cntx,G.p(i).pin(j).cnty],...
            G.p(i).pin(j).a0, G.da, G.n_path, G.r_path);

        
    end

        % debug-plot
%     figure(6); imbo4(G.p(i).im);    hold on
%     for j=1:G.n_pins
%         for k=1:G.n_path
%             plot([G.p(i).pin(j).cntx, G.p(i).pin(j).endpts.x(k)],...
%                  [G.p(i).pin(j).cnty, G.p(i).pin(j).endpts.y(k)],'r')
%         end    
%     end
%     
%%


    % improfile
    mult=1; % multiplied with path length will be # of sampling points
    for j =1:G.n_pins % for 4 centers
        for k=1:G.n_path %for all endpoints
            [G.p(i).pin(j).prof(k,:)]=improfile(G.p(i).im',...
                [G.p(i).pin(j).cntx,G.p(i).pin(j).endpts.x(k)],...
                [G.p(i).pin(j).cnty,G.p(i).pin(j).endpts.y(k)],...
                (G.r_path+1)*mult);

        end
    end
    dist=linspace(0,G.r_path,(G.r_path+1)*mult)/res;
%%
    % deduce thickness
    G.r_wall=5.14; % rod outter wall radius from center
    
    for j=1:4 % centers
        for k=1:G.n_path % angles

            G.p(i).pin(j).thickness(k)=...
                sum(G.p(i).pin(j).prof(k,100:end).*(G.p(i).pin(j).prof(k,100:end)>0.001));
            
            
%             % find walls end % will be for every ray individually later
%             [c,G.p(i).pin(j).idx_wall(k)]=min(abs(dist-G.r_wall));
%             
%             % find the corresponding profile value
%             G.p(i).pin(j).val(k)=G.p(i).pin(j).prof(k,G.p(i).pin(j).idx_wall(k));
%             
%             %find last profile index whos values is 'similar'
%             [c2,G.p(i).pin(j).idx_film(k)]=find(G.p(i).pin(j).prof(k,:)>...
%                 G.p(i).pin(j).val(k),1,'last');
%             G.p(i).pin(j).thickness(k)=sum(G.p(i).pin(j).prof(...
%                 k,G.p(i).pin(j).idx_wall(k):G.p(i).pin(j).idx_film(k)));
           
        end
        G.p(i).pin(j).th_mean=mean(G.p(i).pin(j).thickness);
        
        LFT(iii,j,:)=G.p(i).pin(j).thickness;
        
    end
    iii=iii+1;
    
     

end


    
    
