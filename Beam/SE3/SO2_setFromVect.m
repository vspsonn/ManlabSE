function q = SO2_setFromVect(e1, e2)

s0 = e1(1) + e2(2) + 1.0;
s = -e1(1) - e2(2) + 1.0;

q = zeros(2,1);
if s > s0
    s = 0.5 * sqrt(s + 1.0);
    q(1) = 0.25 / s * (e1(2) - e2(1));
    q(2) = s;
else
    s = 0.5 * sqrt(s0 + 1.0);
    q(1) = s;
    q(2) = 0.25 / s * (e1(2) - e2(1));
end

if q(1) < 0
    q = -q;
end

end