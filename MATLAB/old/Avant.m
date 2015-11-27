function [ sol ] = Avant( m )
format compact; format short e;
clf;
n=64;
tfinal=1;
dx=1/(64+1);
dt=dx;
x=(0:dx:1)';


nombreiteration=ceil(tfinal/dt);
nombreplot=20;  

nombreentreplot=ceil(nombreiteration/nombreplot);
c=1;
u=cos(pi*x);

plot(x,u);ax=axis;axis(ax);hold on;
t=0;
ind=0;

%ratio r
r=c*(dt/(dx^2));


%la diagonale de la matrice A qui doit etre inversé. (Mais il faut pas
%l'inversé...)
princ(1)=1;
princ(2:n+1)=1;
princ(n+2)=1;

%la matrice tri diag
a=diag(princ,0);

i=2;
sol(:,1)=u;
%le membre de gauche
b=ones(n+2,1);

for iteration=1:nombreiteration;

%Ax=b, on crée le vecteur b...
b(1)=u(1)+r*(1*u(2)-1*u(1));

b(2:n+1)=u(2:n+1)+r*(u(3:n+2)-2*u(2:n+1)+u(1:n));

b(n+2)=r*(1*u(n+1)-1*u(n+2))+u(n+2);

u=a\b;
sol(:,i)=u;

i=i+1;
t=t+dt;


ind=ind+1;
if((ind==nombreentreplot)|(iteration==nombreiteration));
ind=0;


un=u/exp(-pi^2*(t));
disp([iteration,t,max(un)]);
plot(x,un);hold on;
xlabel('x');ylabel('u normalise');
titre=sprintf('Temps = %0.5g Iteration = %d',t,iteration);
title(titre);

end;
end

