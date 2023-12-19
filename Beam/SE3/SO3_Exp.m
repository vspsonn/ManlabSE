function q = SO3_Exp(p)

q = zeros(4,1);
q(1) = sqrt(1.0 - 0.25 * (p' * p));
q(2:4) = 0.5 * p;

end