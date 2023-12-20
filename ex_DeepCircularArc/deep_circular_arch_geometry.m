function f = deep_circular_arch_geometry(s, n_dim)

R = 100;
angle_root = 197.5;
angle_tip = -17.5;

angle = angle_root + 0.5 * (s+1) * (angle_tip - angle_root);
ca = cosd(angle);
sa = sind(angle);

if n_dim == 2
    f = zeros(4, 1);
    f(3:4) = [R*ca; R*sa];
    
    tg = [sa; -ca];
    normal_ext = [ca; sa];
    f(1:2) = SO2_setFromVect(tg, normal_ext);
else
    f = zeros(7, 1);
    f(5:7) = [R*ca; R*sa; 0];
    
    tg = [sa; -ca; 0];
    normal_ext = [ca; sa; 0];
    f(1:4) = SO3_setFromVect(tg, normal_ext);    
end

end