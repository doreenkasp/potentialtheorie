function[g]=graviPunkt(x1,x2,y1,y2,z1,dz)

    dx=x2-x1;
    dy=y2-y1;
    z2=z1+dz;
    ds = sqrt(dx^2+dy^2);   
    
    
    r1=sqrt(x1^2+y1^2);
    r2=sqrt(x2^2+y2^2);
    
    C = dy/ds;
    S = dx/ds; 
    
    d1=x1*S+y1*C;
    d2=x2*S+y2*C;
    
 
    P=(x1*y2-x2*y1)/ds;
 
    
    
    R12=sqrt(r1^2+z2^2);
    R11=sqrt(r1^2+z1^2);
    R21=sqrt(r2^2+z1^2);
    R22=sqrt(r2^2+z2^2);
    
    A=acos((x1*x2 + y1*y2)/(r1*r2));
    
    if P>=0
        Sp=1;
    else
        Sp=-1;
    end
 
        
    %Formel 3)
    g = (Sp*A*dz...
       + z2*(atan((z2*d1)/(P*R12)) - atan((z2*d2)/(P*R22)))...
       - z1*(atan((z1*d1)/(P*R11)) - atan((z1*d2)/(P*R21)))...
       - P * log(((R22+d2)/(R12+d1))*((R11+d1)/(R21+d2))));
   
      
   
