%% creates preview- .avi files of existing videos
clear all; 
clc; 
format compact; 
close all;


old_path=pwd;
fpath=('G:\cbolesch\Data\20161208 Campaign 4\');
oldFolder = cd(fpath);

l= dir('*.seq');

imsize=[1024 640];
header_size=2048;

xpix=1:2:imsize(1);
ypix=1:2:imsize(2);

fig=figure(1);
k=0;
for file = l'
    k=k+1;
    disp(file.name)
    fullpath=strcat(file.folder,'\',file.name);
    fid=fopen(fullpath);
    hed = fread(fid,header_size);%the header size migth need to be adjusted depending on image settings
    %to read into one matrix to process fruther with MATLAB comment the above and uncomment this
   % Create a 3D Matrix out of the imported data
    fileSize = file.bytes/2;
    fileSize=floor((fileSize-header_size)/imsize(1)/imsize(2)/2);
    img=uint16((fread(fid,imsize(1)*imsize(2)*fileSize,'uint16')));
    fclose(fid);

    % Create a 3D Matrix out of the imported data
    img=reshape(img,[1024 640 fileSize]);
    
    % reduce memory size by factor 4
    img=single(img(xpix,ypix,:));
    
    % get mean open beam value
    m=single(squeeze(mean(mean(img(25:500,250:315,:))))); 
    m2=mean(m);
    means(k)=m2
    
    a=single(reshape(m/m2,1,1,fileSize));
    b=repmat(a,512,320,1);
    img=img./b;
    
    
    for i=1:size(img,3)
        f_BoCount(i,50,15,5);
        imshow(squeeze(img(:,:,i)),[0,1.3*m2]);
        set(gca,'YDir','normal')    
        drawnow()
        M(i)=getframe(fig);
        
    end
    
    
    viddir='G:\cbolesch\Data\';
    vidname=strcat('2016_12_8__',file.name);

    video=VideoWriter([strcat(viddir,vidname)],'MPEG-4');
    video.FrameRate=30;
    open(video)
    writeVideo(video,M);
    close(video)
    clear img M b

end

cd(oldFolder)
% %%

