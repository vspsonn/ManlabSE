function q = SO3_setFromVect(e1, e2)

e3 = cross(e1, e2);

ind = 0;
s_max = e1(1) + e2(2) + e3(3);

s = e1(1) - e2(2) - e3(3);
if s > s_max
    ind = 1;
    s_max = s;
end

s = -e1(1) + e2(2) - e3(3);
if s > s_max
    ind = 2;
    s_max = s;
end

s = -e1(1) - e2(2) + e3(3);
if s > s_max
    ind = 3;
    s_max = s;
end

s = 0.5 * sqrt(s_max + 1);

q = zeros(4, 1);
if ind == 0
    q(1) = s;
    q(2) = 0.25 / s * (e2(3) - e3(2));
    q(3) = 0.25 / s * (e3(1) - e1(3));
    q(4) = 0.25 / s * (e1(2) - e2(1));
elseif ind == 1
    q(1) = 0.25 / s * (e2(3) - e3(2));
    q(2) = s;
    q(3) = 0.25 / s * (e2(1) + e1(2));
    q(4) = 0.25 / s * (e3(1) + e1(3));
elseif ind == 2
    q(1) = 0.25 / s * (e3(1) - e1(3));
    q(2) = 0.25 / s * (e1(2) + e2(1));
    q(3) = s;
    q(4) = 0.25 / s * (e3(2) + e2(3));
elseif ind == 3
    q(1) = 0.25 / s * (e1(2) - e2(1));
    q(2) = 0.25 / s * (e3(1) + e1(3));
    q(3) = 0.25 / s * (e3(2) + e2(3));
    q(4) = s;
end

if q(1) < 0
    q = -q;
end

end