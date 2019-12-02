function [ output_args ] = f2_FolderCheck(paths)

% takes the folder paths as cell arrays
% if they dont exist, they will be created


for i =1: size(paths,2)
    if exist(paths{i},'dir')~=7
        mkdir(paths{i})
    end
end



; %1

; %2

; %3

; %4

; %5

; %6

; %7

; %8

; %9

; %10;
