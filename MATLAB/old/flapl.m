% etant donne une solution pour le probleme de la membrane rangee
% dans le vecteur phi de taille (nt^2,1) - cette fonction retourne
% le vecteur psi (meme taille que phi) = A phi, avec A la matrice
% de discretisation pour le laplacien sur la membrane.
% note: cette fonction ne forme pas la matrice A !

function psi=flapl(phi);
m=length(phi);
nt=sqrt(m);
n=(nt-1)/2;
h=1/n;
u=reshape(phi,nt,nt);
v=4*u;
% rectangle inferieur 
ic=n+1:nt-1;
jc=2:nt-1;
im=ic-1;
ip=ic+1;
jm=jc-1;
jp=jc+1;
v(ic,jc)=v(ic,jc)-u(im,jc)-u(ip,jc)-u(ic,jm)-u(ic,jp);
% carre en haut a droite
ic=2:n;
jc=n+1:nt-1;
im=ic-1;
ip=ic+1;
jm=jc-1;
jp=jc+1;
v(ic,jc)=v(ic,jc)-u(im,jc)-u(ip,jc)-u(ic,jm)-u(ic,jp);
psi=reshape(v,nt^2,1);
psi=psi/h^2;
