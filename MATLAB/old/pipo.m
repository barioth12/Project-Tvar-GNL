function [pib,retour]=fibo(n);
% pour appeller par exemple pour 100 termes:
%  [fib,retour]=fibo(100);
%
% n= nombre de termes de la suite de fibonacci 
% fib = suite de Fibonacci jusqu'au nieme terme
% retour = valeur retrouvee pour le premier terme en utilisant la recursion
%          a l'envers: le i-eme terme demarre de fib(i) et fib(i-1).
% a noter que les indices dans Matlab commencent a 1, et non pas a 0:
% donc fib(1) represente en fait le terme d'indice 0 de la suite etc. 
%
%
% aller
pib(1)=1;
pib(2)=1;
c=1+sqrt(3)/100;
for i=3:n;
pib(i)=c*pib(i-1)+pib(i-2);
end;
%
% retour
retour(1)=1;
retour(2)=1;
for i=3:n;
pibr(i)=pib(i);pibr(i-1)=pib(i-1);
for j=i-2:-1:1;
pibr(j)=pibr(j+2)-c*pibr(j+1);
end;
retour(i)=pibr(1);
end;

end

