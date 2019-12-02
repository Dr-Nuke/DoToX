function [reg] = f42_RegioFinder(mat,init)
% returns a sorted list of materials corresponding 
% to the 2xn material variable

reg=zeros(1,size(mat,1)+1);
reg(1)=init;
for i = 1:length(mat)
    reg(i+1)=mat(i,mat(i,:)~=reg(i));


end

