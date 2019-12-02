for i=cases
% load file

b=f3_loadSingleVariableMATFile(M.imcname{i});

%extract siogram, add to M
M.sinoslice{i}=squeeze(b(:,45,:));
% clear file
clear b


    
end
