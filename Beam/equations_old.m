function [Rf,dRf] = equations(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]


U = Uf(1:sys.neq+1);      % Main variables
Ua = Uf(sys.neq+2:end);  % Auxiliary variables

q_A = U(1:3);
q_B = U(4:6);
Lambda_BC = U(7:9);
w = U(10);
q0_A = U(11);
q0_B = U(12);
Moment = U(sys.neq+1);

p = Ua(1:3);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue


R(10) = w * w - (1. - 0.25 * p' * p);
R(11) = q0_A * q0_A - (1. - q_A' * q_A);
R(12) = q0_B * q0_B - (1. - q_B' * q_B);


Ra(1:3) = p - (q0_A * q_B - q0_B * q_A - skew(q_A) * q_B);
sigma = sys.parameters.K * p;

loc_gp{1} = [0.];
log_gp{2} = [];

weight_gp{1} = [2.];
weight_gp{2} = [];

n_gp = 1;
ind = 13;
inda = 4;
for k = 1:n_gp
    
    alpha = 0.5 *(loc_gp{n_gp}(k) + 1.);
    
    zk = Ua(inda);
    ck = Ua(inda+1:inda+3);
    bk = U(ind);
    fk = U(ind+1:ind+3);
    
    Ra(inda) = zk - (1. - 0.25 * alpha^2 * p' * p);
    Ra(inda+1:inda+3) = ck - skew(p) * fk;
    
    R(ind) = zk * bk - p' * sigma;
    R(ind+1:ind+3) = zk * fk - (sigma + 0.25 * alpha^2 * bk * p);
    
    ind = ind + 4;
    inda = inda + 4;
    
    R(1:3) = R(1:3) + 0.5 * weight_gp{n_gp}(k) * ( -w * fk - 0.5 * ck );
    R(4:6) = R(4:6) + 0.5 * weight_gp{n_gp}(k) * ( w * fk - 0.5 * ck );
end

R(6) = R(6) + Moment;

R(1:3) = R(1:3) + Lambda_BC;
R(7:9) = q_A;

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
