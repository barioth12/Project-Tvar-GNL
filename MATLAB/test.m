function [ estimator ] = test(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
u=3; sig=4;a=2;b=4;p=4;
%Debute par generer n variable GNL;
n=1000;
X=GNL(u,sig,a,b,p,n);


%t=ones(20,1);
t=0.0001:0.0001:0.002;


estimator=distance(X,t);

end

