function [ sortie ] =GNLcdfconv( x,mu,sigma2,alpha,beta,rho )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


integrandescalaire=@(s)cdf('Normal',x-s,mu*rho,sqrt(sigma2*rho))*GLpdf(s,alpha,beta,rho);

integrande=@(s) arrayfun(integrandescalaire,s);

sortie= integral(integrande,-Inf,Inf,'reltol',1e-8);

end

