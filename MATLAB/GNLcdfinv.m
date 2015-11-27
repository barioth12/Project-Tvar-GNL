function [ sortie ] = GNLcdfinv( x,mu,sigma2,alpha,beta,rho )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
y=x-rho*mu;

cffn= @(s) ((alpha*beta*exp(-sigma2*(s^2)/2))/((alpha-1i*s)*(beta+1i*s)))^rho;

tempfc=@(s) norm(cffn(s))/s*sin(-angle(cffn(s))+s*y);
integrande=@(s)arrayfun(tempfc,s);

sortie=1/2+(1/pi)* integral(integrande,0,Inf,'reltol',1e-8,'AbsTol',1e-8);

end

