function p = SE2_Log(f)

p = zeros(3,1);
p(1) = f(1) * f(3) - f(2) * f(4);
p(2) = f(1) * f(4) + f(2) * f(3);
p(3) = 2.0 * f(2);

end