function [ data ] = datagenerator( p,n,mu,sigma2,alpha,beta,rho )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

min=alphaquantile( p,mu,sigma2,alpha,beta,rho  );

data=zeros(n,1);
compteur=1;
while compteur <=n
    x=GNL(mu,sigma2,alpha,beta,rho,1);
    if x>=min
        data(compteur)=x;
        compteur=compteur+1;
    end
end

end

