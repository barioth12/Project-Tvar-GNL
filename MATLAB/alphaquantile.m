function [ xnew ] = alphaquantile( quantile,mu,sigma2,alpha,beta,rho )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
func=@(s)GNLcdfinv(s,mu,sigma2,alpha,beta,rho)-quantile;
dfunc=@(s) GNLpdfinv(s,mu,sigma2,alpha,beta,rho);

xold=rho*(mu+1/alpha-1/beta);

xnew=xold-func(xold)/dfunc(xold);
compteur=2;
while compteur <100000 && abs(func(xnew))>1e-10
    
    xold=xnew;
    xnew=xold-func(xold)/dfunc(xold);
    
compteur=compteur +1;
end
end

