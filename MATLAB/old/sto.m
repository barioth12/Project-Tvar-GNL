function [ ST ] = sto( t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T=2;    si=1;   r=0.5;  s=2;


dt = .001;

sqrtdt = sqrt(dt);

y=s;
while t<=T
    y = y + dt*(r*y) + si*(cos(2*(T-t))^2)*sqrt(abs(y))*sqrtdt*normrnd(0,1);
    t=t+dt;
end
 ST=y;


end

