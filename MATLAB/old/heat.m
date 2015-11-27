% Ici version inconditionnellement stable avec 
% Crank-Nicholson 
% Experience pour illustrer l'instabilite avec l'equation de
% la chaleur. 
% Je resouds u_t=u_xx avec 0<=t<= tfinal=0.5 et 0<=x<=1
%    condition initiale: u(t=0,x)=sin(pi*x);
%    conditions aux bords: u(t,x=0)=0=u(t,x=1)
%
% la solution exacte est u_ex(t,x)=exp(-pi^2 t) sin(pi x)
%
% je prends ici: dt=dx/2 avec dx=1/(n+1)
%
% je fais regulierement le graphique de u(t,x)/exp(-pi^2 t) qui
% selon la solution exacte devrait etre donne par sin(pi x).
% 

function [sol]= heat(n);
format compact; format short e;
clf;
tfinal=1;
dx=1/(n+1);
p=1/16;

dt=dx/2;
x=(0:dx:1)';


nombreiteration=ceil(tfinal/dt);
nombreplot=20;  

nombreentreplot=ceil(nombreiteration/nombreplot);
c=1+sin(2*pi*(x)/p)/2;
ce=sqrt(3)/2;
u=cos(pi*x);

plot(x,u);ax=axis;axis(ax);hold on;
t=0;
ind=0;
r=(dt/(dx^2));


%la diagonale de la matrice A qui doit etre inversé. (Mais il faut pas
%l'inversé...)
princ(1)=1+r/2*c(1);
princ(2:n+1)=1+r/2*c(2:n+1).*2;
princ(n+2)=1+r/2*c(n+2);

%la diagonale sup.
over(1)=-r/2*c(1);
over(2:n+1)=-r/2*(c(3:n+2)/4-c(1:n)/4+c(2:n+1));
%la diagonal inf
under(1:n)=-r/2*(-c(3:n+2)/4+c(1:n)/4+c(2:n+1));
under(n+1)=-r/2*c(n+2);
%la matrice tri diag
a=diag(princ,0) + diag(under,-1) + diag(over,1);

i=2;
sol(:,1)=u;
%le membre de gauche
b=u;

for iteration=1:nombreiteration;

%Ax=b, on crée le vecteur b...
b(1)=u(1)+r/2*c(1)*(1*u(2)-1*u(1));

b(2:n+1)=u(2:n+1)+r/2*(((c(3:n+2)-c(1:n)).*(u(3:n+2)-u(1:n)))/4+c(2:n+1).*(u(3:n+2)-2*u(2:n+1)+u(1:n)));

b(n+2)=r/2*c(n+2)*(1*u(n+1)-1*u(n+2))+u(n+2);

u=a\b;
sol(:,i)=u;

i=i+1;
t=t+dt;


ind=ind+1;
if((ind==nombreentreplot)|(iteration==nombreiteration));
ind=0;


un=u/(exp(-pi^2*(t)));

disp([iteration,t,max(un)]);
plot(x,un);hold on;
xlabel('x');ylabel('u normalise');
titre=sprintf('Temps = %0.5g Iteration = %d Pour c(x) avec n=%d',t,iteration,n);
title(titre);

end;
end;
