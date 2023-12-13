function [Rf,dRf] = equations_vec(sys,Uf,dUf)
% Equations of the system of the form Rf(Uf) = C + L(Uf) + Q(Uf,Uf).

lambda = Uf(sys.neq+1,:);
u = Uf(1:sys.neq,:);            du = dUf(1:sys.neq,:);
Ua = Uf(sys.neq+2:end,:);       dUa = dUf(sys.neq+2:end,:);

K=sys.parameters.K;

%%% Residues
R =K* u + lambda.*Ua;       dR    = zeros(size(R));
Ra = Ua - exp(u);           dRa = dUa - Ua.*du;

%%% Concatenation
Rf= [R ; Ra];
dRf=[dR;dRa];

end
