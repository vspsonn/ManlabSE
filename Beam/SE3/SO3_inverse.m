function f_inv = SE3_inverse(f)
f_inv = zeros(4,1);
f_inv(1) = f(1);
f_inv(2:4) = f(2:4) * -1;
end