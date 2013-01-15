function [P] = kreis(n,PU,r)

P=zeros(3,n);
gamma=2*pi/n;

for n_i=1:n
   
    P(1,n_i) = r*sin(n_i*gamma)+PU(1);
    P(2,n_i) = r*cos(n_i*gamma)+PU(2);
    P(3,n_i) = PU(3);
    
end
