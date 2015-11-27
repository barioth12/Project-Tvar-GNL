function [ psy ] = fonctionchara( mu,sigma2,alpha,beta,rho,t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
j=1;
for l =t
    psy(j)=((alpha*beta*exp(1i*mu*t(j)-sigma2*t(j)^2/2))/((alpha-1i*t(j))*(beta+1i*t(j))))^rho;
    j=j+1;
end
end

