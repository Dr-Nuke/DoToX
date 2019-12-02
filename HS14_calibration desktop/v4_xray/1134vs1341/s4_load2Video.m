% file info


M.d.fileInfo = dir(M.d.filepath{1});

M.d.fileSize = M.d.fileInfo.bytes;

M.d.nFrames{1}=(M.d.fileSize-M.d.header_size)/...
    M.d.imsize(1)/M.d.imsize(2)/2;

if mod(M.d.nFrames{1},1)~=0
    error('calculated number of video frames is not integer')
end

M.d.fid=fopen(M.d.filepath{1});

% dispose of header
hed = fread(M.d.fid,M.d.header_size);%the header size migth need to be adjusted depending on image settings

%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
 raw1=uint16((fread(M.d.fid,...
    M.d.imsize(1)*M.d.imsize(2)*M.d.fileSize,'uint16'))); %the image
    
fclose(M.d.fid);

% rotate such that the representation is right
raw1=rot90(reshape(raw1,[M.d.imsize(1) M.d.imsize(2) M.d.nFrames{1}]));


save('raw1','raw1')
clear raw1
disp('video 1 done')

%% 2nd file
% file info
M.d.fileInfo = dir(M.d.filepath{2});

M.d.fileSize = M.d.fileInfo.bytes;

M.d.nFrames{2}=(M.d.fileSize-M.d.header_size)/...
    M.d.imsize(1)/M.d.imsize(2)/2;

if mod(M.d.nFrames{2},1)~=0
    error('calculated number of video frames is not integer')
end

M.d.fid=fopen(M.d.filepath{2});

% dispose of header
hed = fread(M.d.fid,M.d.header_size);%the header size migth need to be adjusted depending on image settings

%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
 raw2=uint16((fread(M.d.fid,...
    M.d.imsize(1)*M.d.imsize(2)*M.d.fileSize,'uint16'))); %the image
    
fclose(M.d.fid);

% rotate such that the representation is right
raw2=rot90(reshape(raw2,[M.d.imsize(1) M.d.imsize(2) M.d.nFrames{2}]));


save('raw2','raw2')
disp('video 2 done')
clear raw2
save('M','M')
clear all