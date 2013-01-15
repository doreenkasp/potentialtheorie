clear all
close all
clc

% grad und ordnung festlegen
order=360;   % m max
degree=360;  % n max

% Einlesen
load 'data.dat';

% Konstanten
a=0.6378136300E+07;
r=6378137;
gammaM=3.986005*10^(14);

% Skalierung
delta=200;   % Anzahl Schritte
phi = 0:2*pi/delta:2*pi;
theta = 0:pi/delta:pi;  % theta=0=Nordpol, theta=pi=Suedpol
mu=cos(theta);

% Daten neusortieren
nimpro=data(:,1);
mimpro=data(:,2);
Cnmimpro=data(:,3);
Snmimpro=data(:,4);

% Initialisierung
Snmtmp(1:degree*order+1)=0;
Cnmtmp(1:degree*order+1)=0;
Snm(1:degree+1,1:order+1)=0;
Cnm(1:degree+1,1:order+1)=0;
Ug(1:delta+1,1:delta+1)=0;
Vg(1:delta+1,1:delta+1)=0;
Wg(1:delta+1,1:delta+1)=0;
U(1:delta+1,1:delta+1)=0;
Q(1:delta+1,1:delta+1)=0;

%Laufindex zur Neuordnung
ii=0;

% Bedenke S00=S(1,1) ... S1010 = S(11,11)

% Schleife zur Aussortierung der benutzten Daten
for i=1:length(data(:,1))
    if nimpro(i) < degree+1
        if mimpro(i) < order+1
           ii=ii+1;
           Snmtmp(ii)=Snmimpro(i);
           Cnmtmp(ii)=Cnmimpro(i);
        end
    end
end

% Schleife zur Matrixerstellung (Neusortierung)
n=1:degree+1;
m=1:order+1;
ii=0;
for j=1:order+1
    for i=1:degree+1
        if i >= j
            ii=ii+1;
            Cnm(i,j)=Cnmtmp(ii);
            Snm(i,j)=Snmtmp(ii);
        else
            Cnm(i,j)=0;
            Snm(i,j)=0;
        end
    end
end

% % Kontrollplot Cnm vor Korrektur
% figure(4)
% plot(0:2:10,log(Cnm(1:2:11,1)))


% Korrekturen
Cnm(3,1) = Cnm(3,1) + 0.108262982131 * 10^(-2)/sqrt(5);
Cnm(5,1) = Cnm(5,1) + 0.237091120053 * 10^(-5)/sqrt(9);
Cnm(7,1) = Cnm(7,1) + 0.608346498882 * 10^(-8)/sqrt(13);
Cnm(9,1) = Cnm(9,1) + 0.142681087920 * 10^(-10)/sqrt(17);
Cnm(11,1) = Cnm(11,1) + 0.121439275882 * 10^(-13)/sqrt(21);

% % Kontrollplot Cnm nach Korrektur
% figure(5)
% plot(0:2:10,log(Cnm(1:2:11,1)))

% Vorfaktor
kon = gammaM / r;

k=1;
for n = 0:degree
    P = legendre(n,mu,'sch');
       for j=1:delta+1;
            for m = 0:order
                  if m <= n
                    Q(m+1,:) = P(m+1,:) .* (Cnm(n+1,m+1)*cos(m*phi(j))+Snm(n+1,m+1)*sin(m*phi(j)));
                  end
            end
            Vg(j,:)=sum(Q);
            Wg(j,:)=Wg(j,:)+Vg(j,:);
       end
    U(:,:)=(a/r)^n*Wg(:,:);
    Ug(:,:)=Ug(:,:)+U(:,:);
end
Ug(:,:)=Ug(:,:)*kon;

Ug=Ug';

phi=0:2*pi/delta:2*pi;

figure(1)
contourf(phi/pi*180,-theta/pi*180,Ug,1000,'Linestyle','none')
xlabel('Longitude','Fontsize',24)
ylabel('Latitude','Fontsize',24)
title('Schwerepotential via Legendre Polynome, nmax=10, mmax=10','Fontsize',24)
set(gca,'Fontsize',24)
colormap jet
colorbar

figure(2)
[X,Y,Z] = sphere(length(Ug));
surf(X,Y,Z,Ug)
axis equal
xlabel('\mu','Fontsize',24)
ylabel('\mu','Fontsize',24)
zlabel('\mu','Fontsize',24)
title('Schwerepotential via Legendre Polynome,  nmax=10, mmax=10','Fontsize',24)
set(gca,'Fontsize',24)
colormap jet
colorbar

% % "Beulenplot"
% 
% % Kugelkoordinaten berechnen
% [phi,theta] = meshgrid(phi,theta);
% ord = max(max(abs(Ug)));
% rho = 10^(11) + 2*Ug/ord;
% % Apply spherical coordinate equations
% r = rho.*sin(theta-pi);
% rx = r.*cos(phi);    % spherical coordinate equations
% ry = r.*sin(phi);
% rz = rho.*cos(theta-pi);
% 
% figure(3)
% surf(rx/10^(11),ry/10^(11),rz/10^(11),Ug)
% axis equal
% xlabel('\mu','Fontsize',24)
% ylabel('\mu','Fontsize',24)
% zlabel('\mu','Fontsize',24)
% title('Schwerepotential via Legendre Polynome,  nmax=10, mmax=10','Fontsize',24)
% set(gca,'Fontsize',24)
% colormap jet
% colorbar