function q_3 = SO2_composition(q_1, q_2)

n_cols = size(q_1, 2);

q_3 = zeros(2, n_cols);
q_3(1, :) = q_1(1, :) * q_2(1) - q_1(2, :) * q_2(2);
q_3(2, :) = q_1(2, :) * q_2(1) + q_1(1, :) * q_2(2);

end