function [ sol ] = sigmaimp( S,t,K,T,sigtest,r  )
%Fonction pour calcul� la volatilit� implicite
%sigtest est la vrai valeur de la volatilit� pour obtenir V1, en r�alit� on
%aurais V1 directement.

    %Changement de variable donn� en cours
    u=log((S*exp(r*(T-t)))/K);

    % Point d'inflextion qui sera notre point de depart
    sigma=sqrt((2*u)/(T-t));

    % calcul de V1
    v1=bschole( S,t,K,T,sigtest,r );

    %Changement de variable donn� en cours
    p=sigma*sqrt(T-t);
   
    %Pour sauvegard� les it�rations
    sol=sigma;

    %le calcul de d+
    d1=u/p+p/2;
    
    %la m�thode de newton
    n=1;
    while n<=500;
    %it�ration de la solution
    sol(n)=sigma-(bschole( S,t,K,T,sigma,r )-v1)/(S*sqrt(T-t)*(1/(sqrt(2*pi)))*exp(-(d1^2)/2));
    %test de sortie
    if abs(sol(n)-sigma)<=0.000001
        break;
    end
    sigma=sol(n);
    n=n+1; 
    
    end;
   
    
end

