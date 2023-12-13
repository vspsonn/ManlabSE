function [Rf,dRf] = equations(sys,Uf,dUf)
% Equations of the system of the form Rf(Uf) = C + L(Uf) + Q(Uf,Uf).

% Variables
lambda = Uf(sys.neq+1);
u = Uf(1:sys.neq);      du = dUf(1:sys.neq);
Ua = Uf(sys.neq+2:end); dUa = dUf(sys.neq+2:end);

% Residues
R = zeros(sys.neq,1);                       dR    = zeros(sys.neq,1);
Ra = zeros(sys.neq_aux,1);                  dRa = zeros(sys.neq_aux,1);

% Definition of the auxiliary variables | differentiation of the non-quadratic equations
%  -1   +   v1   -  u1*u2
%   0   +  5u1   +  v1*v2
%   0   +   v3   +          -  tanh(v2) |             dv3 - v4*dv2
%  -1   +   v4   +  v3*v3
%   0   +   v5   +          -  exp (u2) |             dv5 - v5*du2
%   0   +  v6-v5 +  u1*v6
Ra(1) = -1 + Ua(1) - u(1)*u(2);
Ra(2) = 5*u(1) + Ua(1)*Ua(2);
Ra(3) = Ua(3) - tanh(Ua(2));            dRa(3) = dUa(3) - Ua(4)*dUa(2);
Ra(4) = -1 + Ua(4) + Ua(3)*Ua(3);
Ra(5) = Ua(5) - exp(u(2));              dRa(5) = dUa(5) - Ua(5)*du(2);
Ra(6) = Ua(6)-Ua(5) + u(1)*Ua(6);

% Equation of the original system
%   0   +   u1   +  lambda*v6
%   0   +   u2   +  u1*v3
R(1) = u(1) + lambda*Ua(6);
R(2) = u(2) + u(1)*Ua(3);

% Concatenation
Rf=[R ; Ra];
dRf=[dR ; dRa];

end