function [lam,psi]=egmresflap(nt);
format compact;format long e;
n=(nt-1)/2;

% initialise: partout 1 sur le L, 0 au bord et hors du L
phi0=zeros(nt,nt);
phi0(2:nt-1,2:nt-1)=1;
phi0(1:n,1:n)=0;
psi0=reshape(phi0,nt*nt,1);
lami=1;
[max0,ii0]=max(abs(psi0));
%plot
%figure(10);imagesc(phi0);axis('square');pause(1);

%iteration puissance inverse
%===========================
% suggestion :  nombre iterations augmente avec taille
%
itermax=40*(ceil(nt/20))^2;

% voici les iterations
%---------------------
x=psi0;

for iter=1:itermax;

% a completer pour calculer l'inverse de la valeur propre lami et le vecteur propre psi
y=gmresflapl(x);


[lami,p]=max(abs(y));
x=y/lami;

end;

% output final
lam=1/lami;
psi=reshape(psi0,nt,nt);
