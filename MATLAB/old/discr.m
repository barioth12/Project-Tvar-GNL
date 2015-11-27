function [ xi, soli ] = discr( n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h=10/n;
xi=0:h:10;
%solexi=xi.^2+sin(xi)+1;
solexi=asin(tanh(20*(xi/10-1/2))/2);
u=0.2;
% restriction a l'interieur du domaine,
xip=h:h:10-h;
ui=asin(tanh(20*(xip/10-1/2))/2)+u*(1-xip/10);
b=[-1];
for i=1:n-3;
b=[b,-1];
end;
% terme non-homogene dip et membre de droite bip correspondant
dip= -tanh(20*(xip/10-1/2))+exp(ui)-exp(-ui);
for j=1:n-2;
    dip(j)=dip(j)-(exp(ui(j))+exp(-ui(j)))*ui(j);
end;
%dip= tanh(20*(xip/10-1/2))+exp(ui(1))-exp(-ui(1))-(exp(ui(1))+exp(-ui(1)))*ui(1);
bip=-h^2*dip;
% attention aux corrections pour la premiere et la derniere ligne, pour tenir
% compte des conditions aux bords
bip(1)=bip(1)+asin(tanh(20*(0/10-1/2))/2)+0.2;

bip(n-1)= bip(n-1)+asin(tanh(20*(10/10-1/2))/2);
%linearisation ci
ci=(exp(ui)+exp(-ui));
%ci=(exp(ui(1))+exp(-ui(1)));
% construire les elements de la matrice de discretisation
diag0=2+h*h*ci;
diag2=-1-xip(1:n-2).^2-h*xip(1:n-2);
a=diag(diag0)+diag(b,1)+diag(b,-1);
%solution du systeme discretise
solip=a\bip';
soli=[solexi(1)+u,solip',solexi(n+1)];
end