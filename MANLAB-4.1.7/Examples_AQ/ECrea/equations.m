function [Rf,dRf] = equations(sys,Uf,dUf)
% Equations of the system of the form R(U) = C + L(U) + Q(U,U).
% The residu that defines auxiliar variables is computed and then the
% real residu of the system.

u = Uf(1:sys.neq);         % u1,u2
lambda = Uf(sys.neq+1);    % lambda
Ua = Uf(sys.neq+2:end);  % v1, ..., v4

mu = sys.parameters.mu;

% Residues
R = zeros(sys.neq,1);                       dR    = zeros(sys.neq,1);
Ra = zeros(sys.neq_aux,1);                  dRa = zeros(sys.neq_aux,1);

%%% Auxiliary residue
%  ra1 = v1 - 1 - u1 - u1*u1  = 0
%  ra2 = v2 - 1 - u2 - u2*u2  = 0
%  ra3 = 1 - v3*v1            = 0
%  ra4 = 1 - v4*v2            = 0
Ra(1) = Ua(1) -1 - u(1) - u(1)*u(1);
Ra(2) = Ua(2) -1 - u(2) - u(2)*u(2);
Ra(3) = -1 + Ua(1)*Ua(3);
Ra(4) = -1 + Ua(2)*Ua(4);

%%% Main residue
%   r1 = 2 u1 -  u2  + 100 u1 v3 -lambda       = 0
%   r2 = 2 u2 -  u1  + 100 u2 v4 -lambda - mu  = 0
R(1) = 2*u(1) - u(2) + 100*u(1)*Ua(3) - lambda;
R(2) = 2*u(2) - u(1) + 100*u(2)*Ua(4) - lambda -mu;

Rf = [R;Ra];
dRf = [dR;dRa];

end