% solution exemple, version linearisee; version solveur lineaire "lent"
function [xi,soli,solexi]=exemple(n);
%discretisation en x et solution de reference
h=1/n;
xi=0:h:1;
solexi=xi.^2+sin(xi)+1;
% restriction a l'interieur du domaine, 
xip=h:h:1-h;
% terme non-homogene dip et membre de droite bip correspondant
dip=2+6*xip.^2+2*xip.*cos(xip)-(1+xip.^2).*sin(xip);
bip=-h^2*dip;
% attention aux corrections pour la premiere et la derniere ligne, pour tenir
% compte des conditions aux bords
bip(1)=bip(1)+1;
bip(n-1)=bip(n-1)+(2+sin(1))*(1+xip(n-1).^2+h*xip(n-1));
% construire les elements de la matrice de discretisation
diag0=2+2*xip.^2;
diag2=-1-xip(1:n-2).^2-h*xip(1:n-2);
a=diag(diag0)+diag(diag2,1)+diag(diag2,-1);
%solution du systeme discretise
solip=a\bip';
soli=[solexi(1),solip',solexi(n+1)];

