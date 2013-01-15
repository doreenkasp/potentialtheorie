%% POTENZIALTHEORIE -- Aufgabenzettel 2
% Marius Kriegerowski, Moritz Nieschlag, Doreen Kasper und Janina Kammann

clear all
close all
%matlabpool open
%% Deklarationen

drho = 200;
Sm = 1;
gamma = 6.67384E-11;   % Gravitationskonstante
NE = 100;              % Anzahl der Ecken
NS = 100;              % Anzahl der Schichten
PU = [0 0 4];          % Mittelpunkt der Kugel
R = 2;                 % Radius

V = gamma * drho * Sm; 



% Schleife ueber Anzahl der Schichten:
for ns=4:NS
    
    dz = 2*R/(ns-1);             % Schichtdicke

    ne=ns;
    
    % Schleife ueber Anzahl der Ecken
    %parfor ne=4:NE
        
    
        B=kugel(R,ns,ne,PU);    % Kugel erstellen
        gout=0;
        g=0;      
        
        % Schleife ueber Schichten
        for ns_i=1:ns-2
            
            % Schleife ueber Ecken
            for w_i=1:ne
                
                % Koordinaten abgreifen
                x1=B(ns_i,1,w_i);
                y1=B(ns_i,2,w_i);
                z1=B(ns_i,3,w_i);
                x2=B(ns_i,1,mod(w_i,ne)+1);
                y2=B(ns_i,2,mod(w_i,ne)+1);
                
                % Schwere berechnen
                g=graviPunkt(x1,x2,y1,y2,z1,dz);
               
                % Aufsummieren
                gout=gout+g;
            
            end
            
        end
       
      % Abspeichern
      XYZ(ns, ne) = gout*V;
      
      % analytische Loesung berechnen
      g_analyt(ns,ne)=-4/3*pi*gamma*drho*(R-dz)^3*1/((PU(1))^2+(PU(2))^2+PU(3)^2);  
      
      
      gdiff(ns,ne)=(g_analyt(ns,ne)-XYZ(ns,ne))/g_analyt(ns,ne);

  %  end
end

figure(1)
plot(gdiff(50,:))
xlabel('Anzahl der Ecken')
ylabel('Abweichung in Prozent')
title('Prozentuale Abweichung der approximierten Loesung')

figure(2)
plot(gdiff(50,:))
xlabel('Anzahl der Schichten')
ylabel('Abweichung in Prozent')
title('Prozentuale Abweichung der approximierten Loesung')