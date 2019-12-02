function [ phan ] = f4_SquarePhan( l )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here


h=l/2

x=[h,h,-h,-h,h];
y=[-h,h,h,-h,-h];

phan.xy=[x;y];
phan.c=[1,2,3,4];
end

