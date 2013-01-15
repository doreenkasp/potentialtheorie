%clear all
close all

Grad = 3;
Ordnung = 3;

delta = pi/200;
theta = 0 : delta : pi; 
phi = 0 : 2*delta : 2*pi; 
[phi,theta] = meshgrid(phi,theta);

% Legendre ausrechnen
Ymn = legendre(Grad,cos(theta(:,1)),'sch');

%Polynome gesuchter Ordung extrahieren:
Ymn = Ymn(Ordnung+1,:)';

%Grid erstellen:
yy = Ymn;
for kk = 2: size(theta,1)
    yy = [yy Ymn];
end;

%Werte für phi ausrechnen:
yy = yy.*cos(Ordnung*phi);
K=ones(length(yy));
[X,Y,Z]=sphere(length(yy))
surf(X,Y,Z,yy)
