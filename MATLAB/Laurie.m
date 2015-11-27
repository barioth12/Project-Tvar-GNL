function  Laurie(  )
% n=100;
% x=linspace(-3,3,n);
% y=linspace(-3,3,n);
% z=linspace(-3,3,n);
% [X,Y,Z]=ndgrid(x,y,z);
% F=((-(X.^2) .* (Z.^3) -(9/80).*(Y.^2).*(Z.^3)) + ((X.^2) + (9/4).* (Y.^2) + (Z.^2)-1).^3);
% isosurface(F,0)
% lighting phong
% caxis
% axis equal
% colormap('flag');
% view([55 34]);


% Fuzzy heart-shaped Mandelbrot fractal
% Mariano Beguerisse
% October 2010

% iterations
n = 300;
% resolution
N = 200;
% Create grid
x = linspace(-1, 1, N);
y = linspace(-1.4, .6, N);
[X,Y] = meshgrid(x, y);

Z = X + i*Y;
Zn = Z;

figure
for j=1:n
  % Mandelbrot map with random noise
  Zn = -i*(Zn).^2 +  (rand(N,N).^(1/5)).*Z;
  M = abs(Zn);
  ind1 = find(M<2);
  ind2 = find(M>=2);
  M(ind1) = 0;
  M(ind2) = -1;
  m = 1;
  imagesc(x, y, (M/m)) 
  colormap([1 1 1; 1 0 0])
  set(gca,'YDir','normal')
  title(n-j)
  pause(0.01)
end
title('Je taime Laurie ', 'fontsize', 16)
end

