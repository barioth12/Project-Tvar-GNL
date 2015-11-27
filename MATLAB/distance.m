function [ minimum ] = distance(X,t)
%calcul la distance a minimise.
psyn=empcf(X,t);

[muini,sigma2ini,alphaini,betaini,rhoini]=depart(X);
psy= @(mu,sigma2,alpha,beta,rho) ((alpha*beta*exp(mu*1i*t-sigma2*(t.^2)/2))./((alpha-1i*t).*(beta+1i*t))).^rho;

dist=@(x) sum((real(psyn)-real(psy(x(1),x(2),x(3),x(4),x(5)))).^2+(imag(psyn)-imag(psy(x(1),x(2),x(3),x(4),x(5)))).^2);


%dist=@(x) sum((real(psyn)-real(((x(3)*x(4)*exp(x(1)*1i*t-x(2)*(t.^2)/2))./((x(3)-1i*t).*(x(4)+1i*t))).^x(5))).^2+(imag(psyn)-imag(((x(3)*x(4)*exp(x(1)*1i*t-x(2)*(t.^2)/2))./((x(3)-1i*t).*(x(4)+1i*t))).^x(5))).^2);
%solution de l'optimisation
minimum=fminsearchbnd(dist,[muini,sigma2ini,alphaini,betaini,rhoini],[inf,0,0,0,0],[inf,inf,inf,inf,inf]);
%minimum=fminsearch(dist,[2,4,3,4,4]);
%vrai distance
%minimum=dist([2,4,3,4,4]);

%resultat donner par optimisation
%minimum=dist([ 4.0477,    1.2368,    0.9179,   12.1967 ,  -2.5572]);
end

