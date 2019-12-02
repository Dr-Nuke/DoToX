function M = f4_loadVideo(M)
% loads the video file
% M = the main structure

% file info


M.d.fileInfo = dir(M.d.filepath);

M.d.fileSize = M.d.fileInfo.bytes;

M.d.nFrames=floor((M.d.fileSize-M.d.header_size)/...
    M.d.imsize(1)/M.d.imsize(2)/2);

M.d.fid=fopen(M.d.filepath);



% dispose of header
hed = fread(M.d.fid,M.d.header_size);%the header size migth need to be adjusted depending on image settings




%to read into one matrix to process fruther with MATLAB comment the above and uncomment this
 M.raw=uint16((fread(M.d.fid,...
    M.d.imsize(1)*M.d.imsize(2)*M.d.fileSize,'uint16'))); %the image
    
fclose(M.d.fid);

% rotate such that the representation is right
M.raw=rot90(reshape(M.raw,[M.d.imsize(1) M.d.imsize(2) M.d.nFrames]));


end

