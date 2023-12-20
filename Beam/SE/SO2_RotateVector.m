function x = SO2_RotateVector(q, X)

n_cols = size(X, 2);

s = 2. * q(2);
c = 1 - q(2) * s;
s = s * q(1);

x = zeros(2, n_cols);
x(1, :) = c*X(1, :) - s*X(2, :);
x(2, :) = s*X(1, :) + c*X(2, :);

end