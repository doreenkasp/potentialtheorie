clear all;
close all;
clc;

%% Deklaration der Variablen und Konstanten
rad=2;   %Radius der Kugel
drho=200;    %Dichteunterschied Kugel-Umgebungsgestein
G=6.67384*10^-11;   %Gravitationskonstante
%sm=1;    %1 wenn Schwerpunkt Prisma unterhalb Fieldpoint, sonst -1
z_oben=2;
l=5;   %Anzahl der Schichten (ungerade)
delta=21;   %ungerade
fpx=linspace(-5,5,delta);
fpy=linspace(-5,5,delta);

g_analyt=zeros(length(fpx),length(fpx));

for j=1:delta
for i=1:delta
%g_analyt(i,j) = G*drho*4/3*pi*rad^3*(z_oben+rad)/(sqrt(fpx(i)^2+fpy(j)^2+(z_oben+rad)^2))^3;  
g_analyt(i,j)=2*pi*G*drho*(rad^2-(1/3)*(fpx(i)^2+fpy(j)^2));
end
end

diff=zeros(1,98);
for n=3:100
GG=zeros(length(fpx),length(fpx));
%k=1;   %Anzahl der Polygone
%n=100;   %Anzahl der Polygonseiten
k=n;   %Anzahl der Eckpunkten



%% Definieren der Fieldpoints, in den Vektoren fpx und fpy stehen jeweils die x- und
%y-Koordinaten der Fieldspoints


%Berechnung der Eckpunktkoordinaten
A=2*pi/k;
x=zeros(k,1);
y=zeros(k,1);
z=linspace(-2,-6,l+1);
r=zeros(l,1);

phi=linspace(0+(pi/(l-1)/2),pi-(pi/(l-1))/2,l);

for j=1:l
    r(j)=rad*sin(phi(j));  

for i=1:k
    x(i)=r(j)*cos((i-1)*A);
    y(i)=r(j)*sin((i-1)*A);
end     
end

for i=1:delta
    for j=1:delta
        fpX=fpx(i);
        fpY=fpy(j);
        GG(i,j)=gravity(G,drho,n,x,y,fpX,fpY,l,z);       
    end
end

    %diff(n-2)=abs(-GG(i,j)-g_analyt(i,j));
     diff(n-2)=GG((delta+1)/2,(delta+1)/2)-g_analyt((delta+1)/2,(delta+1)/2);
     %Diff(n-2)=sum(diff);
%      Diff(n-2)=diff;
end

%% Plotten der Ergebnisse
figure
contourf(fpx,fpy,GG,100,'LineStyle','none')
title('Schwerepotential einer Kugel')
colorbar;

figure
contourf(fpx,fpy,g_analyt,100,'LineStyle','none')
title('Schwerepotential einer Kugel analytisch')
colorbar;

figure  
x_plot = 3:100;
plot(x_plot,abs(diff),'k')   
xlim([3 100])
