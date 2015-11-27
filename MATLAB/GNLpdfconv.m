function [ sortie ] = GNLpdfconv( x,mu,sigma2,alpha,beta,rho )
normalpdf=@(t) 1/(sqrt(sigma2*rho*2*pi))*exp((-((t) - mu*rho)^2)/(2*sigma2*rho));

integrandescalaire=@(s) normalpdf(x-s)*GLpdf(s,alpha,beta,rho);

integrande=@(s) arrayfun(integrandescalaire,s);

sortie= integral(integrande,-Inf,Inf,'reltol',1e-8,'AbsTol',1e-8);
end

