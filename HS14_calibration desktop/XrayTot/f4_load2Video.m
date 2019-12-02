function M = f4_load2Video(M)
% loads two video files
% M = the main structure

% file info
for i = 1:size(M.d.filepath,2)
M.d.fileInfo{i} = dir(M.d.filepath{i});

M.d.fileSize{i} = M.d.fileInfo{i}.bytes;

M.d.nFrames{i}=floor((M.d.fileSize{i}-M.d.header_size)/...
    M.d.imsize(1)/M.d.imsize(2)/2);
M.d.fid{i}=fopen(M.d.filepath{1});

hed = fread(M.d.fid{1},M.d.header_size);%the header size migth need to be adjusted depending on image settings

end



% dispose of header

%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
M.raw1=uint16((fread(M.d.fid{i},...
    M.d.imsize(1)*M.d.imsize(2)*M.d.fileSize{i},'uint16'))); %the image
    
fclose(M.d.fid{1});

% rotate such that the representation is right
M.raw1=rot90(reshape(M.raw1,[M.d.imsize(1) M.d.imsize(2) M.d.nFrames{1}]));

raw1=M.raw{1};
save('raw1','raw1')
M=rmfields{M,raw{1}};
clear raw1

hed = fread(M.d.fid{2},M.d.header_size);%the header size migth need to be adjusted depending on image settings

%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
 M.raw2=uint16((fread(M.d.fid{i},...
    M.d.imsize(1)*M.d.imsize(2)*M.d.fileSize{2},'uint16'))); %the image
    
fclose(M.d.fid{2});

% rotate such that the representation is right
M.raw2=rot90(reshape(M.raw{2},[M.d.imsize(1) M.d.imsize(2) M.d.nFrames{i}]));

raw2=M.raw2;
save('raw2','raw2')
clear raw2
load('raw1')
M.raw1=raw1;
clear raw1

end

