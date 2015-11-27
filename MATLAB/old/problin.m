function [ xi soli ] = problin( n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h=1/n;
xi=0:h:10;
%solexi=xi.^2+sin(xi)+1;
solexi=asin(tanh(20*(xi/10-1/2))/2);
u=0.2;
% restriction a l'interieur du domaine,
xip=h:h:1-h;
% terme non-homogene dip et membre de droite bip correspondant
lip= tanh(20*(xip/10-1/2));
bip=-h^2*lip;
% attention aux corrections pour la premiere et la derniere ligne, pour tenir
% compte des conditions aux bords
bip(1)=bip(1)+tanh(20*(0/10-1/2))+0.2;

bip(n-1)= bip(n-1)+tanh(20*(10/10-1/2))*(1+xip(n-1).^2+h*xip(n-1));
% construire les elements de la matrice de discretisation
diag0=2+2*xip.^2;
diag2=-1-xip(1:n-2).^2-h*xip(1:n-2);
a=diag(diag0)+diag(diag2,1)+diag(diag2,-1);
%solution du systeme discretise
solip=a\bip';
soli=[solexi(1)+u,solip',solexi(n+1)];
for i=1:m;
    

end;
end

