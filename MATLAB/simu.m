function [ int ] = simu(p, mu,sigma2,alpha,beta,rho)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=5000;
data=datagenerator( p,n,mu,sigma2,alpha,beta,rho );
int=bootstrap( data,n ); 
end

