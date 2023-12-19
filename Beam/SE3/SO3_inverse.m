function q_inv = SO3_inverse(q)

q_inv = q;
q_inv(2:4) = q_inv(2:4) * -1;

end