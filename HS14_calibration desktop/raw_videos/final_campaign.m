%% creates preview- .avi files of existing videos
clear all; 
clc; 
format compact; 
close all;

%needed for bocount
addpath('C:/Users/cbolesch/Desktop/HS14_calibration desktop')



cd('C:\Users\cbolesch\Desktop\HS14_calibration desktop\raw_videos')
old_path=pwd;
pathbase='G:\cbolesch\Data\201712xx final campaign\';
paths={'20171212 final campaign\',...
        '20171213 final campaign\',...
        '20171215 final campaign\',...
        '20171218 final campaign\'};

header_size=2048; % vid file header size in bytes    
    
for i=3:4%1:size(paths,2)
    dpath=strcat(pathbase,paths{i}); %directory path
    oldFolder = cd(dpath);
    VidConv(dpath,0,header_size);
    
end
cd('C:\Users\cbolesch\Desktop\HS14_calibration desktop\raw_videos')
 
function VidConv(dpath,hzflag,header_size) 
% converts a xray video into a small size preview video
% dpath = folder path
    
    l= dir('*.seq');
        for file = l'

         
            imsize=[1024 640];
            dose=[460,635,120,1010];     % x start value
 
            if contains(file.name,'85') | hzflag==1
                imsize=imsize/2;
                dose=dose/2;
            end
            xpix=1:2:imsize(1);
            ypix=1:2:imsize(2);
            
            [img,fileSize]=VidRead(file,imsize,header_size);
            % Create a 3D Matrix out of the imported data
            img=reshape(img,[imsize(1) imsize(2) fileSize]);

%% manual fixes
%             if all(file.name=='20171213_annularFlow_30Hz.seq')
%                 img(:,:,1:40)=[];
%                 fileSize=size(img,3);
%             end
            
            
            
            % get mean open beam value
            dosey=dose(1):dose(2);
            dosex=dose(3):dose(4);
            m=single(squeeze(mean(mean(img(dosex,dosey,:))))); 
            m2=mean(m);
            
            a=single(reshape(m/m2,1,1,fileSize));
            b=repmat(a,imsize(1),imsize(2),1);
            img=single(img)./b;
                   
            % reduce memory size by factor 4
            img=single(img(xpix,ypix,:));
            
            % crop one pixel line 
            img=img(4:end-2,2:end,:);

            % manually flip the image
            img=flip(img,1);
            
            wfile=strcat(file.folder,'\','conv_',file.name);
            vmin=min(img(:));
            vmax=max(img(:));
            VidMaker(img,wfile,vmin,vmax)
        end
    a=1;
        
end


function [img,fileSize]=VidRead(file,imsize,header_size)
    disp(file.name)
    fullpath=strcat(file.folder,'\',file.name);
    fid=fopen(fullpath);
    hed = fread(fid,header_size);%the header size migth need to be adjusted depending on image settings
    %to read into one matrix to process fruther with MATLAB comment the above and uncomment this
    % Create a 3D Matrix out of the imported data
    fileSize = file.bytes;
    fileSize=(fileSize-header_size)/imsize(1)/imsize(2)/2;
    if mod(fileSize,1)~=0
        error(strcat(fullpath,': file size is not integer'))
    end

    img=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16=>uint16')));
    fclose(fid);
end

function VidMaker(img,wfile,vmin,vmax)
%% makes a video of the x*y*z stack img
% img = x*y*z array
    fig=figure(1);
    for i=1:size(img,3)
        f_BoCount(i,50,15,5);
        imshow(squeeze(img(:,:,i)),[vmin,vmax]);
        %set(gca,'YDir','normal')    
        drawnow()
        M(i)=getframe(fig);

    end

    video=VideoWriter(wfile,'MPEG-4');
    video.FrameRate=30;
    open(video)
    writeVideo(video,M);
    close(video)

end