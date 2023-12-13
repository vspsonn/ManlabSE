function f = SO3_Exp(p)
f = zeros(4,1);

f(1) = sqrt(1.0 - 0.25 * (p' * p));
f(2:4) = 0.5 * p;

end