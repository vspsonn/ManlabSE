function [Rf,dRa,Forcing] = equations(sys,t,Uf,dUf,d2Uf)
% Equations of the system of the form Rf(Uf) = C + L(Uf) + Q(Uf,Uf).

u = Uf(1:sys.nz);         % x,y
du = dUf(1:sys.nz);       % x',y'
Ua = Uf(sys.nz+1:end);    % v,r
lambda = Uf(end);         % lambda

R = zeros(sys.nz,1);
Ra = zeros(sys.nz_aux,1);
dRa = zeros(sys.nz_aux,1);

alpha = sys.parameters.alpha;

% Equations of the system
% x' = y
% y' = -alpha y + v*r - x
R(1) = u(2) - du(1);
R(2) = -alpha*u(2) + Ua(1)*Ua(2) - u(1) - du(2);

% Auxiliary variables
%     0 = v - lambda y
%     0 = r - 1 + x^2
Ra(1) = Ua(1) -lambda*u(2);
Ra(2) = Ua(2) -1 - u(1)^2 ;

%%% Concatenation
Rf = [R;Ra];


%%% Forcing
Forcing = zeros(length(t),sys.nz_tot);

end