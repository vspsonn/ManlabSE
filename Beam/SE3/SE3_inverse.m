function f_inv = SE3_inverse(f)

f_inv = zeros(7,1);
f_inv(1) = f(1);
f_inv(2:4) = f(2:4) * -1;
f_inv(5:7) = - (eye(3) - 2.0 * f(1) * skew(f(2:4)) + 2.0 * skew(f(2:4)) * skew(f(2:4))) * f(5:7);

end