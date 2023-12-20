function f = SE2_Exp(p)

f = zeros(4,1);
f(1) = sqrt(1.0 - 0.25 * p(3) * p(3));
f(2) = 0.5 * p(3);
sca = 1/f(1) * (1.0 - f(2)^2);
f(3) = sca * p(1) - f(2) * p(2);
f(4) = sca * p(2) + f(2) * p(1);

end