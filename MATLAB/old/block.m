function [ norm ] = block( A)
%algo 2.1 de l'article du tp4

%déclare certaine variable qui nous serons utiles

%n est la taille de la matrix
n=size(A);
n=n(1);
%le nombre d'itération maximum.
s=n/5;

ai=eye(n);

format rat;

%vecteur des colones non visité.
nonvis=1:n;
i=0;
x=mean(ai(:,nonvis)')';
%repeat
while i<s
% Déclaration de mon vecteur x de norme 1 ayant seulement des valeur sur
% les position non visité.
%x=mean(ai(:,nonvis)')';

%Le vecteur y défini tel que dans l'algo
y=A*x;

%le vecteur sign défini tel que dans l'algo
si=sign(y);

%le vecteur z defini tel que dans l'algi
z=A'*si;

%une condition de sortie si on trouve la norme exact.
if max(abs(z))<= (z'*x)
       norm=sum(abs(y));              
       return;
end

%trouve la valeur maximal de z et ca position
[m2,j2]=max(abs(z(nonvis)));
%change le vecteur nonvis pour retiré la position du z_max trouvé
jnew=nonvis(j2);
nonvis=nonvis(nonvis~=jnew);
x=ai(:,jnew);
%pour sortir de la boucle.
i=i+1;
end


% heuristic choise intended to "pick out" any large element of A in those
% case where such elements fail to be revealed during the algorithm.
i=1;
while i<=n
    b(i)=((-1)^(i+1))*(1+(i-1)/(n-1));
    i=i+1;
end

% Pose la norm comme etant le max(||y|| , ||Ab|| / ||b||  )
 norm=max(sum(abs(y)),(sum(abs(A*b'))/(sum(abs((b'))))));
end