function  f2_parforSave(path,var)
% this function only does save  a variable; 
% it is neccessary because the parfor will not allow
% 'save'; however a separate function containing save works.

save(path,'var');

end

