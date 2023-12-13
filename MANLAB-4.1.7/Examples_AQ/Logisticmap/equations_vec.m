function [Rf] = equations_vec(sys,Uf)
% Quadratic system of the form R(U) = C + L(U) + Q(U,U) + Lh(dU) + Qh(U,dU)

mu = Uf(sys.neq+1,:);
muinv = Uf(sys.neq+2,:);

x = Uf(1:sys.neq,:);

R = muinv.*x(2:end,:) - x(1:end-1,:).*(1-x(1:end-1,:));

Ra = mu.*muinv - 1;

Rf= [x(1,:) - sys.parameters.x0 ; R ; Ra];

end
