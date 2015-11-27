function [ sci ] = bootstrap( data,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

capable=@(x) mean(x);
sci = bootci(n,{capable,data},'alpha',0.01,'type','student') ;
end

