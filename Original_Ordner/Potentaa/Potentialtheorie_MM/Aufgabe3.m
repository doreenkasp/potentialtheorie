%% POTENZIALTHEORIE
% Marius Kriegerowski, Moritz Nieschlag, Doreen Kasper und Janina Kammann

clear all
close all
%matlabpool open
%% Deklarationen


drho = 200;
Sm = 1;
gamma = 6.67384E-11;   % Gravitationskonstante
ne = 100;               % Anzahl der Ecken
ns = 100;               % Anzahl der Schichten
PU = [0 0 4];          % Mittelpunkt der Kugel
R = 2;                 % Radius

Xrange=-15:0.5:15;
Yrange=Xrange;

k=1;

%for ne=5:5:200

    ns=ne;
    dz = 2*R/(ns-1);       % Schichtdicke

    [X,Y] = meshgrid(Xrange, Yrange);
    V = gamma * drho * Sm; 

    XYZ = zeros(length(Xrange), length(Yrange));
	
	% Erzeugung der Kugel:
    B=kugel(R,ns,ne,PU);

for x_i=1:length(Xrange)
    asd=x_i
    for y_i=1:length(Yrange)

        gout=0;
        g=0;
        
        parfor ns_i=1:ns-2
            
            for w_i=1:ne

            
            x1=B(ns_i,1,w_i)-Xrange(x_i);
            y1=B(ns_i,2,w_i)-Yrange(y_i);
            z1=B(ns_i,3,w_i);
            
            x2=B(ns_i,1,mod(w_i,ne)+1)-Xrange(x_i);
            y2=B(ns_i,2,mod(w_i,ne)+1)-Yrange(y_i);
            %z2=B(ns_i+1,3,mod(w_i,ne)+1);
            
            g=graviPunkt(x1,x2,y1,y2,z1,dz);

            gout=gout+g;
            
            end
            
        end
        
       
       XYZ(x_i, y_i) = gout*V;
      g_analyt(x_i,y_i)=-4/3*pi*gamma*drho*(R-dz)^3*1/...
          ((PU(1)-Xrange(x_i))^2+(PU(2)-Yrange(y_i))^2+PU(3)^2);
   
      gdiff(x_i,y_i)=g_analyt(x_i,y_i)-XYZ(x_i, y_i);
      
    end
end

diff(k,1)=gdiff(30,30);
diff(k,2)=ns;

k=k+1

%matlabpool close

figure(1)
contourf(XYZ)
xlabel('X-Achse [m]')
ylabel('Y-Achse [m]')
title('Schwerefeld [m/s^2] einer Kugel, 100-Schichten, 100-eckigen Polygonen')
colorbar
axis('square')

figure(2)
contourf(g_analyt)
xlabel('X-Achse [m]')
ylabel('Y-Achse [m]')
title('Schwerefeld [m/s^2] einer Kugel, analytische Loesung')
colorbar
axis('square')

figure(3)
plot(gdiff(30,:))
xlabel('X-Achse [m]')
ylabel('Differenz')
%ylabel('Y-Achse [m]')
title('Differenz [m/s^2] zwischen genaeherter und analytischer Loesung')
colorbar
axis('square')

figure(4)
surf(XYZ)
xlabel('X-Achse [m]')
ylabel('Y-Achse [m]')
title('Schwerefeld [m/s^2] einer Kugel, 100-Schichten, 100-eckigen Polygonen')
colorbar
axis([-30 30 -30 30 -3E-8 0])


figure(5)
surf(g_analyt)
xlabel('X-Achse [m]')
ylabel('Y-Achse [m]')
title('Schwerefeld [m/s^2] einer Kugel, analytische Loesung')
colorbar
axis('square')

figure(6)
plot(diff(:,2),diff(:,1))
xlabel('Anzahl der Schichten/Ecken')
ylabel('Differenz Analytische und Gerechnete Loesung')
title('Fehlerentwicklung in Abhaengigkeit von der Anzahl der Ecken/Schichten')