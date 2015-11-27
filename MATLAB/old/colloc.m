function [phic,x]=colloc(n,m);
% base sur GO 6.4, collocation
% attention: n defini comme dans le livre: x(1)=0; x(n)=L, h=L/(n-1); 
% avec n points de collocation (points aux bords + points interieurs) 
% de x(1) a x(n);

L=10;
u0=0.2;
h=L/(n-1);
x=0:h:L;x=x';

% je resouds probleme linearise autour de phi_lin=phid+ u0 *(1-x/L)
% je resouds par collocation pour fonction psi
% avec phic=psi + droite qui passe par (x0,phi0) et (xL,phiL)
% et psi(0)=psi(L)=0
%
% on a donc: psi''- 2 cosh(phi_lin) psi = membre de droite 
% membre de droite = 2 sinh(phi_lin) - 2 cosh(phi_lin) (phi_lin-droite) -c(x)

% membre de droite

cx=tanh(20*x/L-10);
phid=asinh(cx/2);
phi_lin=phid+u0*(1-x/L);

%On refait les solutions m fois dans la méthode de newton pour avoir la
%forme non linéaire
phic=0;
j=1;
while j<=m 
cx=tanh(20*x/L-10);
if j>1
phi_lin=phic;
end

droite=phi_lin(1)*(1-x/L)+phi_lin(n)*x/L;
fx=2*sinh(phi_lin)-2*cosh(phi_lin).*(phi_lin-droite)-cx;
qx=-2*cosh(phi_lin);

% elles correspondent aux formules du livre, corrigees !

princ=zeros(1,n);
princ(1)=-9/h^2;
princ(2)=-27/2/h^2+qx(2)*15/4;
princ(n-1)=-27/2/h^2+qx(n-1)*15/4;
princ(n)=-9/h^2;
for i=3:n-2;princ(i)=-3/h^2+qx(i);end;

under=zeros(1,n-1);
under(1)=3/2/h^2+qx(2)/4;
under(2)=6/h^2+qx(3);
under(n-1)=9/h^2;
for i=3:n-2;under(i)=3/2/h^2+qx(i)/4;end;

over=zeros(1,n-1);
over(1)=9/h^2;
over(n-2)=6/h^2+qx(n-2); % over
over(n-1)=3/2/h^2+qx(n-1)/4;
for i=2:n-3;over(i)=3/2/h^2+qx(i)/4;end;


a=diag(princ,0) + diag(under,-1) + diag(over,1);

% resolution 
ci=a\fx;

% evaluation de la fonction aux noeuds
psi(1,1)=0;
psi(2,1)=ci(1)/4+ci(2)*(4-1/4)+ci(3)/4;
psi(3,1)=ci(2)  +ci(3)        +ci(4)/4;

psi(n-2,1)=ci(n-1)+ci(n-2)    +ci(n-3)/4;
psi(n-1,1)=ci(n)/4+ci(n-1)*(4-1/4)+ci(n-2)/4;
psi(n,1)=0;

psi(4:n-3,1)=ci(3:n-4,1)/4+ci(4:n-3,1)+ci(5:n-2,1)/4;

phic=psi+phi_lin(1)*(1-x/L)+phi_lin(n)*x/L;

j=j+1;
end

end