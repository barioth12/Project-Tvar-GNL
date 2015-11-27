function [  ] = GenerationDeDataProjet(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nombreIteration = 20;

compteur =0;
n=5000 ;

fileID = fopen('data3.txt','w');
for i=0:nombreIteration
 i
r=rand(1,5);
data=datagenerator(0.999,n,r(1),r(2),r(3),r(4),r(5));
boot=bootstrap(data,n);
tvarexp=tvar(0.999,r(1),r(2),r(3),r(4),r(5));

if(tvarexp>=boot(1)&&tvarexp<=boot(2))
    compteur=compteur+1;
    in='oui'
else
    in='non'
end   

var=alphaquantile(0.999,r(1),r(2),r(3),r(4),r(5));

fprintf(fileID,'mu=%2i,sigma^2=%2i,alpha=%2i,beta=%2i,rho=%2i$ Z var=%3.4f Z tvar=%3.4f Z P(%3.4f,%3.4f)P %s\r\n',r(1),r(2),r(3),r(4),r(5),var,tvarexp,boot(1),boot(2),in);
%%fprintf(fileID,'\mu=%2i,\sigma^2=%2i,\alpha=%2i,\beta=%2i,\rho=%2i$ Z %3.4f Z %3.4f Z P(%3.4f,%3.4f)P\\ \hline %s\r\n',r(1),r(2),r(3),r(4),r(5),var,tvarexp,boot(1),boot(2),in);

end
fprintf(fileID,'Compteur=%10i et NombreIteration=%10i',compteur,nombreIteration);
fclose(fileID);
