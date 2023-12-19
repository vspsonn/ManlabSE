function f_inv = SE2_inverse(f)

f_inv = zeros(4,1);
f_inv(1) = f(1);
f_inv(2) = f(2) * -1;
sca_0 = (1. - 2. * f(2)^2);
sca_1 = - 2. * f(1) * f(2);
f_inv(3) = -(sca_0 * f(3) - sca_1 * f(4));
f_inv(4) = -(sca_0 * f(4) + sca_1 * f(3));

end