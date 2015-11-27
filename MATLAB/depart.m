function [ mu,sigma2,alpha,beta,rho ] = depart( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

k=1:5;
m=k;
for j=1:5
    m(j)=mean((X).^j);
end
k(2)=m(2)-m(1)^2;
k(3)=2*m(1)^3-3*m(1)*m(2)+m(3);
k(4)=-6*m(1)^4+12*m(1)^2*m(2)-3*m(2)^2-4*m(1)*m(3)+m(4);
k(5)=24*m(1)^5-60*m(1)^3*m(2)+20*m(1)^2*m(3)-10*m(2)*m(3)+5*m(1)*(6*m(2)^2-m(4))+m(5);
    
       F=@(y) [12*k(3)*(y(1)^(-5)-y(2)^(-5))-k(5)*(y(1)^(-3)-y(2)^(-3));
            3*k(3)*(y(1)^(-4)+y(2)^(-4))-k(4)*(y(1)^(-3)-y(2)^(-3))];
    

memoire = fsolve(F,1:2);
alpha=(memoire(1));
beta=(memoire(2));

rho=k(3)/(alpha^-4+alpha^-4)/6;
sigma2=k(2)/rho-alpha^-2-beta^-2;
mu=k(1)/rho-alpha^-1+beta^-1;
end

