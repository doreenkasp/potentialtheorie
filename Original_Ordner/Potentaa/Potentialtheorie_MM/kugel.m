function[A]=kugel(R, ns, ne, PM)

dz =2*R/(ns-1);

for ns_i=1:ns-1
    
    PU(1)=PM(1);
    PU(2)=PM(2);
    PU(3)=PM(3)+R-ns_i*dz;
    
    r=sqrt(R^2-(R-ns_i*dz)^2);

    A(ns_i,:,:) = kreis(ne,PU,r); 
end