clear all 
close all

ns=100
ne=100

B=kugel(2,ns,ne,[0 0 4]);

figure
hold on

for dim=1:ns-2
    
    for w_i=1:ne
        w(w_i)=B(dim,1,w_i);
        v(w_i)=B(dim,2,w_i);
        u(w_i)=B(dim,3,w_i);
    end

    plot3(w,v,u,'marker','o')

end

axis('square')
grid on