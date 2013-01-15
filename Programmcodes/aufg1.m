clear all
close all
%% m<=n
n=9; %Grad 
m=0; %Ordnung

delta = pi/60;
theta = 0 : delta : pi;
phi = 0 : 2*delta : 2*pi;
[phi, theta] = meshgrid(phi,theta);

Pmn = legendre(n, cos(theta(:,1)), 'sch');
Pmn = Pmn(m+1,:)';

for i = 1: size(phi,1)
    Pmnphi(:,i) = Pmn(:);
end;

F = Pmnphi.*(cos(m*phi)+sin(m*phi));

figure
[x,y,z]=sphere(length(F));
surf(x,y,z,F)
axis square
xlabel('x')
ylabel('y')
zlabel('z')
title('Zonales Legendre Polynom mit n=9 und m=0')
colorbar;
%% m<=n
n=6; %Grad 
m=6; %Ordnung

delta = pi/60;
theta = 0 : delta : pi;
phi = 0 : 2*delta : 2*pi;
[phi, theta] = meshgrid(phi,theta);

Pmn = legendre(n, cos(theta(:,1)), 'sch');
Pmn = Pmn(m+1,:)';

for i = 1: size(phi,1)
    Pmnphi(:,i) = Pmn(:);
end;

F = Pmnphi.*(cos(m*phi)+sin(m*phi));

figure
[x,y,z]=sphere(length(F));
surf(x,y,z,F)
axis square
xlabel('x')
ylabel('y')
zlabel('z')
title('Sektorales Legendre Polynom mit n=m=6')
colorbar;

%% m<=n
n=9; %Grad 
m=3; %Ordnung

delta = pi/60;
theta = 0 : delta : pi;
phi = 0 : 2*delta : 2*pi;
[phi, theta] = meshgrid(phi,theta);

Pmn = legendre(n, cos(theta(:,1)), 'sch');
Pmn = Pmn(m+1,:)';

for i = 1: size(phi,1)
    Pmnphi(:,i) = Pmn(:);
end;

F = Pmnphi.*(cos(m*phi)+sin(m*phi));

figure
[x,y,z]=sphere(length(F));
surf(x,y,z,F)
axis square
xlabel('x')
ylabel('y')
zlabel('z')
title('Tesserales Legendre Polynom mit n=9 und m=3')
colorbar;