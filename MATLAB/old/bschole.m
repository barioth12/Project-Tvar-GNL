function [ sol ] = bschole( S,t,K,T,sigma,r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
    u=log((S*exp(r*(T-t)))/K);
    p=sigma*sqrt(T-t);
    d1=u/p+p/2;
    d2=u/p-p/2;
    sol=S*(normcdf(d1)-exp(-u)*normcdf(d2));
    
    %d1=(log(S/K)+(r+((sigma^2)/2))*(T-t))/(sigma*sqrt(T-t));
    %d2=(log(S/K)+(r-((sigma^2)/2))*(T-t))/(sigma*sqrt(T-t));
    %sol=(S*normcdf(d1))-((K*exp(-r*(T-t)))*normcdf(d2));
    

end

