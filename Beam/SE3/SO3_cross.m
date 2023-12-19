function X = SO3_cross(x, y)

n_cols = size(x, 2);
X = zeros(3, n_cols);
X(1, :) = x(2, :) .* y(3, :) - x(3, :) .* y(2, :);
X(2, :) =  x(3, :) .* y(1, :) - x(1, :) .* y(3, :);
X(3, :) = x(1, :) .* y(2, :) - x(2, :) .* y(1, :);

end