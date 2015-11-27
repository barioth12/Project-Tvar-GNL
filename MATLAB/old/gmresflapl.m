% solution pour psi du systeme:  A psi = phi
% pour un vecteur phi donne, avec A
% la matrice de discretisation correspondant
% au laplacien sur la membrane.
% phi et psi sont des vecteurs de tailles (nt^2,1)
% Note: cette fonction resoud le systeme sans former
% la matrice A, en utilisant une methode iterative.
function psi=gmresflapl(phi);
[psi,flag]=gmres(@flapl,phi,[],1.e-12);
