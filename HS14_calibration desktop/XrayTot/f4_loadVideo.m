function M = f4_loadVideo(M,i)
    M.d.FilePath{i}=strcat(M.d.rawfolder,M.d.filenames{i});
    M.d.FileInfo{i} = dir(M.d.FilePath{i});
    M.d.FileSize(i) = M.d.FileInfo{i}.bytes;
    M.d.nFrames(i)=(M.d.FileSize(i)-M.d.header_size)/...
        M.d.imsize(1)/M.d.imsize(2)/2;
    if mod(M.d.nFrames(i),1) ~= 0
        error('file size is not coherent with pixel sizes and header');    
    end
     
    % open the file
    M.d.fid(i)=fopen(M.d.FilePath{i});

    % dispose of header
    hed = fread(M.d.fid(i),M.d.header_size);%the header size migth need to be adjusted depending on image settings

    %to read into one matrix to process fruther with MATLAB comment the above and uncomment this
    M.raw=uint16((fread(M.d.fid(i),...
        M.d.imsize(1)*M.d.imsize(2)*M.d.nFrames(i),'uint16'))); %the image
    fclose(M.d.fid(i));

    % rotate such that the representation is right
    M.raw=rot90(reshape(M.raw,[M.d.imsize(1) M.d.imsize(2) M.d.nFrames(i)]));


end

