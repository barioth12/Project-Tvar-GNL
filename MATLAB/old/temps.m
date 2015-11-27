function [ time ] = temps( n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
time=0;
q=1;
while q< 2
matri=(randn(n,n));
[l,u,p]=lu(matri);

shittime=cputime;
b=blockinv(l,u,p);
timer=cputime-shittime;
time=time+timer;
q=q+1;
end;



end
% 
% n=6000;
% i=1;
% x=6000:1000:18000;
% while n<=18000
% matri=(randn(n,n));
% [l,u,p]=lu(matri);
% time=cputime;
% b=blockinv(l,u);
% t(i)=cputime-time;
% n=n+1000;
% i=i+1;
% end


% n=200;
% i=1;
% x=200:25:800;
% while n<=800
% matri=(randn(n,n));
% 
% time=cputime;
% b=condest(matri);
% t(i)=cputime-time;
% n=n+25;
% i=i+1;
% end
