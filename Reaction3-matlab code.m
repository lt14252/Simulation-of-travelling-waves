function [OUT1,OUT2,OUT3,OUT4] = reaction(V1,V2,V3,V4,t)
    initial = [V1,V2,V3,V4];
    
    k1 = 1800e-3;%mM-1s-1, k1 reaction rate: H2O2+Amplexred---Resofurin;
    k2 = 500e-3;%mM-1s-1,  k2 reaction rate: H2O2+Resofurin--Resazurin;
    sol = ode23s(@func,[0 t],initial);
    OUT1 = sol.y(1,end); 
    OUT2 = sol.y(2,end);
    OUT3 = sol.y(3,end);
    OUT4 = sol.y(4,end); 
 
%     plot(sol.y,t);
%     legend('A','B','C')
    function dydt = func(t,y)
        dydt = zeros(4,1);
        dydt(1) = -k1*y(1)*y(2)-k2*y(1)*y(3);
        dydt(2) = -k1*y(1)*y(2);
        dydt(3) = k1*y(1)*y(2)-k2*y(1)*y(3);
        dydt(4) = k2*y(1)*y(3);
    end
end