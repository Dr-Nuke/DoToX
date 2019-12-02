function [a] = f_normalize(a)
a=(a-min(a(:)))/(max(a(:))-min(a(:)));

end