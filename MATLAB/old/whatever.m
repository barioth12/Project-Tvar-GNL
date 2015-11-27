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

function heat(n);
format compact; format short e;
clf;
tfinal=1;
dx=1/(n+1);
dt=dx/2;
x=(0:dx:1)';
p=1/n;

nombreiteration=ceil(tfinal/dt);
nombreplot=20;
nombreentreplot=ceil(nombreiteration/nombreplot);
c=1+sin(2*pi*(x+dx/2)*p)/2;
u=sin(pi*x);
plot(x,u);ax=axis;axis(ax);hold on;
t=0;
ind=0;
e=ones(n+2,1);
a=spdiags([-dt/dx^2/2*e (dt/dx^2+1)*e -dt/dx^2/2*e],-1:1,n+2,n+2);
a(1,2)=0;a(n+2,n+1)=0; a(1,1)=1; a(n+2,n+2)=1;
for iteration=1:nombreiteration;
t=t+dt;
b=u;
b(2:n+1)=u(2:n+1)+dt/dx^2/2*(c(3:n+2)*(u(3:n+2)-u(2:n+1))-c(2:n+1)*(u(2:n+1)+u(1:n)));
u=a\b;
ind=ind+1;
if((ind==nombreentreplot)|(iteration==nombreiteration));
ind=0;
un=u/exp(-pi^2*t);
disp([iteration,t,max(un)]);
plot(x,u/exp(-pi^2*t));hold on;
xlabel('x');ylabel('u normalise');
titre=sprintf('Temps = %0.5g Iteration = %d',t,iteration);
title(titre);
pause;
end;
end;
