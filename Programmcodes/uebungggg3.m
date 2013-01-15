%% POTENZIALTHEORIE -- Aufgabenzettel 2
% Marius Kriegerowski, Moritz Nieschlag, Doreen Kasper und Janina Kammann

clear all
close all

%% Deklarationen

drho = 200;
Sm = 1;
gamma = 6.67384E-11;   % Gravitationskonstante
NE = 50;              % Anzahl der Ecken
NS = 50;              % Anzahl der Schichten
PU = [0 0 4];         % Mittelpunkt der Kugel
R = 2;                 % Radius

Xrange=0:NE;
Yrange=0:NS;
[X,Y] = meshgrid(Xrange, Yrange);

for ns=4:NS
 EnEs=ns
    for ne=4:ns

dz = 2*R/(ns-1);             % Schichtdicke

V = gamma * drho * Sm; 

B=kugel(R,ns,ne,PU);

        
        gout=0;
        g=0;
        parfor ns_i=1:ns-2
            
            for w_i=1:ne
                
            x1=B(ns_i,1,w_i);
            y1=B(ns_i,2,w_i);
            z1=B(ns_i,3,w_i);
            
            x2=B(ns_i,1,mod(w_i,ne)+1);
            y2=B(ns_i,2,mod(w_i,ne)+1);
       %     z2=B(ns_i,3,mod(w_i,ne)+1);
            
            g=graviPunkt(x1,x2,y1,y2,z1,dz);
%     if g ~= NaN
%         wurstx1 = x1
%         wurstx2 = x2
%         wy1 = y1
%         wy2 = y2
%         wz1 = z1
%  
%     end
            gout=gout+g;
            
            end
            
        end

            
       
      XYZ(ns,ne) = gout*V;
      g_analyt(ns,ne)=-4/3*pi*gamma*drho*(R-dz)^3*1/(PU(1)^2+PU(2)^2+PU(3)^2);
   
      %gdiff(ns,ne) = 100*(XYZ(ns,ne)-g_analyt(ns,ne))/g_analyt(ns,ne);
        gdiff(ns,ne) = (XYZ(ns,ne)-g_analyt(ns,ne))/g_analyt(ns,ne);
      
      
    end
end

figure
plot(gdiff(50,:))
title('ns=50 Abhaengigkeit von Anz. d Ecken')

figure
plot(gdiff(:,50))
title('ne=50 Abhaengigkeit von Anz. d Schichten')
