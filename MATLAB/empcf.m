function [ psyn ] = empcf( X,t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n=length(t);
psyn=zeros(1,n);
for j=1:n
    psyn(1,j)=mean(exp(1i*t(j)*X));
end
end

