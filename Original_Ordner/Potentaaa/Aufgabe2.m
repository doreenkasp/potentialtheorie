clear all
close all

%% Einlesen des Datenfiles 
filename = 'osu89a-mod.gfc';
[A,delimiterOut]=importdata(filename)

gamma = 6.67384E-11;% Gravitationskonstante
M = 5.974E24;       % Erdmasse
r = 6371000.;       % mittlerer Erdradius
%a = 149.6E6;        % grosse Halbachse
V = gamma*M/r;      % Vorfaktor
a = 0.6378136460E+07;

delta = pi/100;
theta = 0 : delta : pi;     % Breite
phi = 0 : 2*delta : 2*pi;   % Länge
%Anlegen der Koeffizienenmatrizen i  denen jedem Grad7Ordnung Koeffizienten
%zuegeordnet werden
for line_i=1:length(A.data)
    Cnm(A.data(line_i,1)+1,A.data(line_i,2)+1)=A.data(line_i,3);
    Snm(A.data(line_i,1)+1,A.data(line_i,2)+1)=A.data(line_i,4);
end

%Korrekturen
Cnm(3,1) = Cnm(3,1) + 0.108262982131 * 10^(-2)/sqrt(5);
Cnm(5,1) = Cnm(5,1) - 0.237091120053 * 10^(-5)/sqrt(9);
Cnm(7,1) = Cnm(7,1) + 0.608346498882 * 10^(-8)/sqrt(13);
Cnm(9,1) = Cnm(9,1) - 0.142681087920 * 10^(-10)/sqrt(17);
Cnm(11,1) = Cnm(11,1) + 0.121439275882 * 10^(-13)/sqrt(21);

n=200; %maximale Anzahl der Koeffizienten
yy=zeros(length(theta));
xy=zeros(length(theta));
asdf(1:length(phi),1:length(phi))=0;
for n_i=0:n %Schleife über n
    %Berechnung der Legendre Polynome (siehe Aufg.1)
    test=n_i
    Pnm = legendre(n_i,cos(theta),'sch'); 
    for phi_i=1:length(phi)
       for m_i=0:n_i
          %Berechnung des Potenzials nach Formel 2.1
            tmp(m_i+1,:) = Pnm(m_i+1,:).*(Cnm(n_i+1,m_i+1)*...
                cos(m_i*phi(phi_i))+Snm(n_i+1,m_i+1)*sin(m_i*phi(phi_i))); 
       end
       %Anlegen der Summen
       asdf(phi_i,:)=sum(tmp);
       xy(phi_i,:)=xy(phi_i,:)+asdf(phi_i,:);
    end
    qwer(:,:) = xy(:,:)*((a/r)^n_i);
    yy(:,:)=yy(:,:)+qwer(:,:);
end
yy(:,:)=yy(:,:)*V;
yy=yy';

%Normierung der Colourbar
delta_value=min(min(yy));
 %U=U./delta_value-1; Größemordnung der Abweic
 yy=yy-delta_value;
 delta_value
 
 
%Graphische Darstellungen
[phiplot,thetaplot]=meshgrid(phi,theta);
figure(1)
%contourf(theta,phi, yy, 100, 'linestyle','None')
contourf(phiplot*180/pi,thetaplot*180/pi, yy, 100,'linestyle','None')
set(gca,'YDir','reverse')
colorbar
%axissquare
title('Schwerepotential der Erde mit nmax=200')
xlabel('Phi');
ylabel('Theta');



















% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=100;
% yy=zeros(length(theta));
% xy=zeros(length(theta));
% asdf(1:length(phi),1:length(phi))=0;
% for n_i=0:n
%     Pnm = legendre(n_i,cos(theta),'sch');
%     for phi_i=1:length(phi)
%        for m_i=0:n_i
%           
%             tmp(m_i+1,:) = Pnm(m_i+1,:).*(Cnm(n_i+1,m_i+1)*cos(m_i*phi(phi_i))+...
%                      Snm(n_i+1,m_i+1)*sin(m_i*phi(phi_i))); 
%        end
%        asdf(phi_i,:)=sum(tmp);
%        xy(phi_i,:)=xy(phi_i,:)+asdf(phi_i,:);
%     end
%     qwer(:,:) = xy(:,:)*((a/r)^n_i);
%     yy(:,:)=yy(:,:)+qwer(:,:);
% end
% yy(:,:)=yy(:,:)*V;
% yy=yy';