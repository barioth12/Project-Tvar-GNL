function [ X ] = GNL( mu,sigma2,alpha,beta,rho,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X=rho*mu+sqrt(sigma2*rho)*normrnd(0,1,n,1)+(1/alpha)*gamrnd(rho,1,n,1) -(1/beta)*gamrnd(rho,1,n,1);

end

