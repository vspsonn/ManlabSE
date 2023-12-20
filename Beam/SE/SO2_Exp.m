function q = SO2_Exp(p)

q = zeros(2,1);
q(1) = sqrt(1.0 - 0.25 * p * p);
q(2) = 0.5 * p;

end