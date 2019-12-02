function [table_out] = f3_xls(path)
% reads and modifies the xls tables 
table_out=readtable(path);
names={'process','ProjNr','ImNr','exp_start','exp_end',...
    'Dose','good','X','sliceNr','Angle','fileNr'};
table_out.Properties.VariableNames=names;

addnames={'ImType','OldFileName','NewFileName'};
addcol=zeros(size(table_out,1),1);
cellcol=cell(size(table_out,1),1);
table_out=[table_out,table(cellcol,cellcol,cellcol,'VariableNames',addnames)];

%
end








%

