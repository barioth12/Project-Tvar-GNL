function [ y ] = pdfcurve(x,mu,sigma2,alpha,beta,rho )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y=x;
compteur=1;
for i=x
    y(compteur)=GNLpdfconv(i,mu,sigma2,alpha,beta,rho);
    compteur=compteur+1;
end

end

