function [ norm,j ] = blockinv( l,u )
%algorithme 2.1 de l'article donner pour le tp4 adapté pour trouver la
%norme de l'inverse d'une matrice.

n=size(l);
n=n(1);
ai=eye(n);
format rat;
nonvis=1:n;
j=0;
x=mean(ai(:,nonvis)')';

i=1;
while i<=n
  b(i)=((-1)^(i+1))*(1+(i-1)/(n-1));
 i=i+1;
end

fore=l\b';
sorti=(sum(abs((u\fore))/(sum(abs((b'))))));

%repeat

while j<=1


 %y=A\x;
fore=l\x;
y=u\fore;

si=sign(y);

%z=A'\si
fore=u'\si;
z=l'\fore;

if max(abs(z))<= (z'*x) &&j>0
       j=j+1;
       %norm=max(sum(abs(y)),sorti);
       norm=sum(abs(y));
       return;
        
end
[m2,j2]=max(abs(z(nonvis)));
jnew=nonvis(j2);
nonvis=nonvis(nonvis~=jnew);
x=ai(:,jnew);
j=j+1;
end


norm=sum(abs(y));
%Pose la norm comme etant le max(||y|| , ||Ab|| / ||b||  )
norm=max(sum(abs(y)),sorti);
end
% 
% i=0;
% compteur1=zeros(3,1);
% compteur10=zeros(3,1);
% compteur50=zeros(3,1);
% sumi=zeros(3,1);
% temps=zeros(3,1);
% while i< 1000
% ai=eye(100);
% format rat;
% matri=(randn(100,100));
% [l,u,p]=lu(matri);
% c=norm(matri,1);
% invmatri=inv(p*matri);
% a=norm(invmatri,1);
% shittime=cputime;
% b=blockinv(l,u);
% temps(2)=temps(2)+cputime-shittime;
% %methode 2
% sumi(2)=sumi(2)+abs(b*c)/(a*c);
% if min(b/a,a/b)>=0.99
% compteur1(2)=compteur1(2)+1;
% end
% if min(b/a,a/b)>=0.9
% compteur10(2)=compteur10(2)+1;
% end
% if min(b/a,a/b)>=0.5
% compteur50(2)=compteur50(2)+1;
% end
% %methode 3
% shittime=cputime;
% condest(matri,1);
% temps(3)=temps(3)-shittime+cputime;
% sumi(3)=sumi(3)+abs(condest(l*u,1))/(a*c);
%  if min((c*a)/condest(matri,4),condest(matri,1)/(c*a))>=0.99
% compteur1(3)=compteur1(3)+1;
% end
%  if min((c*a)/condest(matri,4),condest(matri,1)/(c*a))>=0.9
% compteur10(3)=compteur10(3)+1;
% end
% if min((c*a)/condest(matri,4),condest(matri,1)/(c*a))>=0.5
% compteur50(3)=compteur50(3)+1;
% end
% %methode 1
% j=1;
% sol=0;
% shittime=cputime;
% while j<=20
% val = sum(abs(invmatri(:,j)));
% if val >= sol
% sol = val;
% end
% j=j+1;
% end
% temps(1)=temps(1)-shittime+cputime;
% sumi(1)=sumi(1)+abs(sol*c)/(a*c);
% if min(sol/a,a/sol)>=0.99
% compteur1(1)=compteur1(1)+1;
% end
% if min(sol/a,a/sol)>=0.9
% compteur10(1)=compteur10(1)+1;
% end
% if min(sol/a,a/sol)>=0.5
% compteur50(1)=compteur50(1)+1;
% end
% 
% 
% i=i+1;
% end
% 
% i=1;
% while i<= 1000
% matri=randn(100,100);
% [l,u,p]=lu(matri);
% [a,c]=blockinv(l,u);
% b(i)=c;
% i=i+1;
% end
