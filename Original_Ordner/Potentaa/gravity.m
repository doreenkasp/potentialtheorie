function[g]=gravity(G,drho,n,x,y,fpX,fpY,l,z)

g=zeros(l,1);
g1=zeros(l,1);
for j=1:l   %Schleife über die Schichten
    
    z1=z(j);   %Oberkante Polygon
    z2=z(j+1);   %Unterkante Polygon

for i=1:n      %Schleife über die Polygonkanten
            if i>1
            x1=x(i)-fpX;
            x2=x(i-1)-fpX;
            y1=y(i)-fpY;
            y2=y(i-1)-fpY;
            else
            x1=x(n)-fpX;
            x2=x(i)-fpX;
            y1=y(n)-fpY;
            y2=y(i)-fpY;  
            end
        
r1=sqrt(x1^2+y1^2);
r2=sqrt(x2^2+y2^2);
ds=sqrt((x1-x2)^2+(y1-y2)^2);
A=acos((x1*x2+y1*y2)/(r1*r2));
P=(x1*y2-x2*y1)/(ds);

    if P<0
        sp=-1;
    else 
        sp=1;
    end

d1=x1*(x2-x1)/ds+y1*(y2-y1)/ds;
d2=x2*(x2-x1)/ds+y2*(y2-y1)/ds;

R12=sqrt(r1^2+z2^2);
R22=sqrt(r2^2+z2^2);
R11=sqrt(r1^2+z1^2);
R21=sqrt(r2^2+z1^2);

g1(l)=sp*A*(z2-z1)+z2*(atan((z2*d1)/(P*R12))-atan((z2*d2)/(P*R22)))-z1*(atan((z1*d1)/(P*R11))-atan((z1*d2)/(P*R21)))-P*log((R22+d2)/(R12+d1)*(R11+d1)/(R21+d2));    
g(l)=g(l)+g1(l);
end

g(l)=g(l)*G*drho;
end
g=sum(g);
