function [ var_out ] = s2_boLoad(path)
% loads the single variable of a mat file
a=load(path);
var_out=a.var;

end

