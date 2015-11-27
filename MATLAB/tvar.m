function [ value ] = tvar( p,mu,sigma2,alpha,beta,rho )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


pip=alphaquantile(p,mu,sigma2,alpha,beta,rho);

myfun=@(x) x*GNLpdfconv(x,mu,sigma2,alpha,beta,rho);


integrande=@(s)arrayfun(myfun,s);

t=alphaquantile(1,mu,sigma2,alpha,beta,rho);
x=920:1000;
%plot(x,integrande(x))

%value1= integral(integrande,pip,t,'AbsTol',1e-8,'reltol',1e-8)/(1-GNLcdfconv(pip,mu,sigma2,alpha,beta,rho));
%value2= quadgk(integrande,t,Inf,'AbsTol',1e-8,'RelTol',1e-8)/(1-GNLcdfconv(pip,mu,sigma2,alpha,beta,rho));
%value=value1+value2;
value=integral(integrande,pip,Inf,'AbsTol',1e-8,'reltol',1e-8)/(1-GNLcdfconv(pip,mu,sigma2,alpha,beta,rho));
end

