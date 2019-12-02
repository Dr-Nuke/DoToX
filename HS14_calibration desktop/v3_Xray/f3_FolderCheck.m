function [] = f3_FolderCheck(path )
%F3_FOLDERCHECK Summary of this function goes here
%   Detailed explanation goes here

if exist(path,'dir')~=7
    mkdir(path)
    disp(['created ',path])
else
   disp(['already existing ',path]) 
end
end

