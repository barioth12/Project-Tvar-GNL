function [phig,x]=galerkin(n,m);
% base sur GO 6.4, methode de Galerkin
% attention: (n+2) points en tout, n points interieurs (comme dans le livre)

L=10;
u0=0.2;
h=L/(n+1);
x=0:h:L;x=x';


% je resouds probleme linearise autour de phi_lin=phid+ u0 *(1-x/L)
% je resouds au sens de galerkin pour fonction psi
% avec phig=psi + droite qui passe par (x0,phi0) et (xL,phiL)
% et psi(0)=psi(L)=0
%
% dans cette version - je suppose que je connais phi_lin uniquement
% aux noeuds (aux sommets des fonctions "chapeau"). cette version
% devrait etre plus simple a transformer en methode de Newton
%
% on a donc: psi''- 2 cosh(phi_lin) psi = membre de droite 
% membre de droite = 2 sinh(phi_lin) - 2 cosh(phi_lin) (phi_lin-droite) -c(x)

% conditions aux bords;
cx0=tanh(-10);
cxL=tanh(10);
phid0=asinh(cx0/2);
phidL=asinh(cxL/2);
phi0=phid0+u0;
phiL=phidL;

% construction de la solution de depart aux sommets des chapeaux.
% pour iterer, substituer a phi_lin la derniere solution calculee
cx=tanh(20*x/L-10);
phid=asinh(cx/2);
phi_lin=phid+u0*(1-x/L);

% integration aux 2 abscisses de Gauss x_inti
% dans chaque intervalle, celle a gauche = x_intim et celle a droite = x_intip
% evaluation des fonctions a integrer en ces points la.
% comme la solution phi_lin autour de laquelle je linearise est connue
% uniquement en x, pas aux points de Gauss - j'interpole lineairement par
% morceaux pour trouver les fonctions en ces points la

%On trouve les 3 points dans l'intervale (x(i-1);x(i))
x_intim=x(1:n+1)+h/2*(1-sqrt(3/5));
x_intip=x(1:n+1)+h/2*(1+sqrt(3/5));
x_mid=x(1:n+1)+h/2;

%Calcul la valeur de la fonction c(x) au 3 point dans
%l'intervale[x(i-1),x(i)]

cxm=tanh(20*x_intim/L-10);
cxp=tanh(20*x_intip/L-10);
cxmid=tanh(20*x_mid/L-10);

%Les poids pour l'interpolation
poids1=1/2*(1-sqrt(3/5));
poids2=1/2*(1+sqrt(3/5));

%Les valeur de l'interpolation pour phi(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
phi_linm=poids1*phi_lin(1:n+1)+poids2*phi_lin(2:n+2);
phi_linp=poids2*phi_lin(1:n+1)+poids1*phi_lin(2:n+2);
phi_mid=1/2*phi_lin(1:n+1)+1/2*phi_lin(2:n+2);

%Les valeur de l'interpolation pour droite(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
droitem=phi0*(1-x_intim/L)+phiL*x_intim/L;
droitep=phi0*(1-x_intip/L)+phiL*x_intip/L;
droite_mid=phi0*(1-x_mid/L)+phiL*x_mid/L;

%Les valeur de l'interpolation pour f(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
fxm=2*sinh(phi_linm)-2*cosh(phi_linm).*(phi_linm-droitem)-cxm;
fxp=2*sinh(phi_linp)-2*cosh(phi_linp).*(phi_linp-droitep)-cxp;
fx_mid=2*sinh(phi_mid)-2*cosh(phi_mid).*(phi_mid-droite_mid)-cxmid;

%Les valeur de l'interpolation pour q(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
qxm=-2*cosh(phi_linm);
qxp=-2*cosh(phi_linp);
qxmid=-2*cosh(phi_mid);

% calcul des integrales Ri,Qi,Si,fi1,fi2
%w1 et w3 les "weights" Wi de la méthode, mais ici w2=w1.

w1=5/9;
w3=8/9;

rim=qxm(2:n+1).*(x_intim(2:n+1)-x(3:n+2)).^2;
rip=qxp(2:n+1).*(x_intip(2:n+1)-x(3:n+2)).^2;
rimid=qxmid(2:n+1).*(x_mid(2:n+1)-x(3:n+2)).^2;
ri=h/2*((rim+rip)+rimid);

qim=qxm(1:n).*(x_intim(1:n)-x(1:n)).^2;
qip=qxp(1:n).*(x_intip(1:n)-x(1:n)).^2;
qimid=qxmid(1:n).*(x_mid(1:n)-x(1:n)).^2;
qi=h/2*(w1*(qim+qip)+w3*qimid);

sim(2:n)=qxm(2:n).*(x_intim(2:n)-x(2:n)).*(x_intim(2:n)-x(3:n+1));
sip(2:n)=qxp(2:n).*(x_intip(2:n)-x(2:n)).*(x_intip(2:n)-x(3:n+1));
simid(2:n)=qxmid(2:n).*(x_mid(2:n)-x(2:n)).*(x_mid(2:n)-x(3:n+1));
si(2:n)=h/2*(w1*sim(2:n)+w1*sip(2:n)+w3*simid(2:n));

fi1m=fxm(1:n).*(x_intim(1:n)-x(1:n));
fi1p=fxp(1:n).*(x_intip(1:n)-x(1:n));
fi1mid=fx_mid(1:n).*(x_mid(1:n)-x(1:n));

fi2m=-fxm(2:n+1).*(x_intim(2:n+1)-x(3:n+2));
fi2p=-fxp(2:n+1).*(x_intip(2:n+1)-x(3:n+2));
fi2mid=fx_mid(2:n+1).*(x_mid(2:n+1)-x(3:n+2));

fi=1/2*((fi1p+fi1m+fi2p+fi2m)+(fi1mid+fi2mid));

%construction non-optimale de la matrice a
%construction non-optimale de la matrice a
princ=zeros(1,n);
under=zeros(1,n-1);
over=zeros(1,n-1);
for i=1:n; princ(i)=(-2*h+qi(i)+ri(i))/h^2;end;
for i=1:n-1; over(i)=(h-si(i+1))/h^2;end;
for i=1:n-1; under(i)=(h-si(i))/h^2;end;

a=diag(princ,0) + diag(under,-1) + diag(over,1);
% resolution (a optimiser pour matrice tridiag)

% resolution (a optimiser pour matrice tridiag)
ci=a\fi;

% evaluation de la fonction aux noeuds
psi(1,1)=0;
psi(2:n+1,1)=ci;
psi(n+2,1)=0;

phig=psi+phi0*(1-x/L)+phiL*x/L;

j=1;
while j<=m
    
% construction de la solution de depart aux sommets des chapeaux.
% pour iterer, substituer a phi_lin la derniere solution calculee
cx=tanh(20*x/L-10);
phid=asinh(cx/2);
phi_lin=phig;

% integration aux 2 abscisses de Gauss x_inti
% dans chaque intervalle, celle a gauche = x_intim et celle a droite = x_intip
% evaluation des fonctions a integrer en ces points la.
% comme la solution phi_lin autour de laquelle je linearise est connue
% uniquement en x, pas aux points de Gauss - j'interpole lineairement par
% morceaux pour trouver les fonctions en ces points la

%On trouve les 3 points dans l'intervale (x(i-1);x(i))
x_intim=x(1:n+1)+h/2*(1-sqrt(3/5));
x_intip=x(1:n+1)+h/2*(1+sqrt(3/5));
x_mid=x(1:n+1)+h/2;

%Calcul la valeur de la fonction c(x) au 3 point dans
%l'intervale[x(i-1),x(i)]

cxm=tanh(20*x_intim/L-10);
cxp=tanh(20*x_intip/L-10);
cxmid=tanh(20*x_mid/L-10);

%Les poids pour l'interpolation
poids1=1/2*(1-sqrt(3/5));
poids2=1/2*(1+sqrt(3/5));

%Les valeur de l'interpolation pour phi(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
phi_linm=poids1*phi_lin(1:n+1)+poids2*phi_lin(2:n+2);
phi_linp=poids2*phi_lin(1:n+1)+poids1*phi_lin(2:n+2);
phi_mid=1/2*phi_lin(1:n+1)+1/2*phi_lin(2:n+2);

%Les valeur de l'interpolation pour droite(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
droitem=phi0*(1-x_intim/L)+phiL*x_intim/L;
droitep=phi0*(1-x_intip/L)+phiL*x_intip/L;
droite_mid=phi0*(1-x_mid/L)+phiL*x_mid/L;

%Les valeur de l'interpolation pour f(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
fxm=2*sinh(phi_linm)-2*cosh(phi_linm).*(phi_linm-droitem)-cxm;
fxp=2*sinh(phi_linp)-2*cosh(phi_linp).*(phi_linp-droitep)-cxp;
fx_mid=2*sinh(phi_mid)-2*cosh(phi_mid).*(phi_mid-droite_mid)-cxmid;

%Les valeur de l'interpolation pour q(x) avec x les 3 points dans
%l'intervale[x(i-1),x(i)]
qxm=-2*cosh(phi_linm);
qxp=-2*cosh(phi_linp);
qxmid=-2*cosh(phi_mid);

% calcul des integrales Ri,Qi,Si,fi1,fi2
%w1 et w3 les "weights" Wi de la méthode, mais ici w2=w1.

w1=5/9;
w3=8/9;

rim=qxm(2:n+1).*(x_intim(2:n+1)-x(3:n+2)).^2;
rip=qxp(2:n+1).*(x_intip(2:n+1)-x(3:n+2)).^2;
rimid=qxmid(2:n+1).*(x_mid(2:n+1)-x(3:n+2)).^2;
ri=h/2*(w1*(rim+rip)+w3*rimid);

qim=qxm(1:n).*(x_intim(1:n)-x(1:n)).^2;
qip=qxp(1:n).*(x_intip(1:n)-x(1:n)).^2;
qimid=qxmid(1:n).*(x_mid(1:n)-x(1:n)).^2;
qi=h/2*(w1*(qim+qip)+w3*qimid);

sim(2:n)=qxm(2:n).*(x_intim(2:n)-x(2:n)).*(x_intim(2:n)-x(3:n+1));
sip(2:n)=qxp(2:n).*(x_intip(2:n)-x(2:n)).*(x_intip(2:n)-x(3:n+1));
simid(2:n)=qxmid(2:n).*(x_mid(2:n)-x(2:n)).*(x_mid(2:n)-x(3:n+1));
si(2:n)=h/2*(w1*sim(2:n)+w1*sip(2:n)+w3*simid(2:n));

fi1m=fxm(1:n).*(x_intim(1:n)-x(1:n));
fi1p=fxp(1:n).*(x_intip(1:n)-x(1:n));
fi1mid=fx_mid(1:n).*(x_mid(1:n)-x(1:n));

fi2m=-fxm(2:n+1).*(x_intim(2:n+1)-x(3:n+2));
fi2p=-fxp(2:n+1).*(x_intip(2:n+1)-x(3:n+2));
fi2mid=fx_mid(2:n+1).*(x_mid(2:n+1)-x(3:n+2));

fi=1/2*((fi1p+fi1m+fi2p+fi2m)+(fi1mid+fi2mid));

%construction non-optimale de la matrice a
%construction non-optimale de la matrice a
princ=zeros(1,n);
under=zeros(1,n-1);
over=zeros(1,n-1);
for i=1:n; princ(i)=(-2*h+qi(i)+ri(i))/h^2;end;
for i=1:n-1; over(i)=(h-si(i+1))/h^2;end;
for i=1:n-1; under(i)=(h-si(i))/h^2;end;

a=diag(princ,0) + diag(under,-1) + diag(over,1);
% resolution (a optimiser pour matrice tridiag)

% resolution (a optimiser pour matrice tridiag)
ci=a\fi;

% evaluation de la fonction aux noeuds
psi(1,1)=0;
psi(2:n+1,1)=ci;
psi(n+2,1)=0;

phig=psi+phi0*(1-x/L)+phiL*x/L;
j=j+1;
end
end




% Interpolation de la solution phi(x(i-1),x(i)) pour les 3 point dans
% lintervale
%
% poids1=1/2*(1-sqrt(3/5));
% poids2=1/2*(1+sqrt(3/5));
% 
% %phi_linm=poids1*phi_lin(1:n+1)+poids2*phi_lin(2:n+2);
% %phi_linp=poids2*phi_lin(1:n+1)+poids1*phi_lin(2:n+2);
% %phi_mid=1/2*phi_lin(1:n+1)+1/2*phi_lin(2:n+2);

