function [ sol ] = sigmaimp( S,t,K,T,sigtest,r  )
%Fonction pour calculé la volatilité implicite
%sigtest est la vrai valeur de la volatilité pour obtenir V1, en réalité on
%aurais V1 directement.

    %Changement de variable donné en cours
    u=log((S*exp(r*(T-t)))/K);

    % Point d'inflextion qui sera notre point de depart
    sigma=sqrt((2*u)/(T-t));

    % calcul de V1
    v1=bschole( S,t,K,T,sigtest,r );

    %Changement de variable donné en cours
    p=sigma*sqrt(T-t);
   
    %Pour sauvegardé les itérations
    sol=sigma;

    %le calcul de d+
    d1=u/p+p/2;
    
    %la méthode de newton
    n=1;
    while n<=500;
    %itération de la solution
    sol(n)=sigma-(bschole( S,t,K,T,sigma,r )-v1)/(S*sqrt(T-t)*(1/(sqrt(2*pi)))*exp(-(d1^2)/2));
    %test de sortie
    if abs(sol(n)-sigma)<=0.000001
        break;
    end
    sigma=sol(n);
    n=n+1; 
    
    end;
   
    
end

