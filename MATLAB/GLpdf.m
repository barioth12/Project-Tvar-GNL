function [ sortie ] = GLpdf( x,alpha,beta,tau )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
lbox1=(tau/2+1/4)*log(alpha*beta)+(beta-alpha)*x/2-log(gamma(tau))-log(pi)/2 -(tau/2-1/4)*log(2);

 box1=exp(lbox1);
 
 lbox2=(tau-1/2)*(log(2*alpha*beta)/2 +log(abs(x))-log(alpha+beta)) ;
 
 lbox3=log(besselk(abs(tau-1/2),abs(x)*(alpha+beta)/2,1))-abs(x)*(alpha+beta)/2;
 
if tau>0.5
    box23=(1/2)*gamma(tau - 1/2)*(4*sqrt(2*alpha*beta)*(alpha + beta)^(-2))^(tau - 1/2);
else
    box23=gamma(0.5 - tau)*2^(-2*tau)*(sqrt(2*alpha*beta))^(tau - 0.5)*(abs(x))^(2*tau - 1);
end
box0=box1*box23;

if x==0
    sortie=box0;
else
    sortie=exp(lbox1+lbox2+lbox3);
end
end

