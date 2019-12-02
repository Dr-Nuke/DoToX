function f2_FolderCheck(paths)

% takes the folder paths as cell arrays
% if they dont exist, they will be created



if exist(paths{i},'dir')~=1
    mkdir(paths{i})
    disp(['created ',paths])
else
   disp(['already existing ',paths]) 
end
end




