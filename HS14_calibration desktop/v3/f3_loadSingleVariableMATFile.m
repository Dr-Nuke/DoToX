function y = loadSingleVariableMATFile(filename)
fprintf('loadind %s .....',filename)
foo = load(filename);
whichVariables = fieldnames(foo);
if numel(whichVariables) == 1
    y = foo.(whichVariables{1});
    fprintf('done\n')
else
    error('the file loaded contains another number of variables than 1. nothing was loaded.') 
end
