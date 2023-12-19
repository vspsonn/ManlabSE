function f = bend45_geometry(s, n_dim)

L = 100;

angle = 0.5 * (s+1) * 0.25 * pi;
ca = cos(angle);
sa = sin(angle);

if n_dim == 2
    f = zeros(4, 1);
    f(3:4) = [L*sa; L*(1-ca)];
    
    tg = [ca; sa];
    normal_ext = [-sa; ca];
    f(1:2) = SO2_setFromVect(tg, normal_ext);
else
    f = zeros(7, 1);
    f(5:7) = [L*sa; L*(1-ca); 0];
    
    tg = [ca; sa; 0];
    normal_ext = [-sa; ca; 0];
    f(1:4) = SO3_setFromVect(tg, normal_ext);    
end

end