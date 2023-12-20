function p = SE3_Log(f)

p = zeros(6,1);
p(4:6) = 2.0 * f(2:4);
p(1:3) = f(1) * f(5:7) + skew(f(5:7)) * f(2:4);

end