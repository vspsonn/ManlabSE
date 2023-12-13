function [Rf,dRf] = equations_rot(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]


U = Uf(1:sys.neq+1);      % Main variables
Ua = Uf(sys.neq+2:end);  % Auxiliary variables

q0_A = U(1);
q_A = U(2:4);
q0_B = U(5);
q_B = U(6:8);

Lambda_BC = U(9:11);
Moment = U(sys.neq+1);

p = Ua(1:3);
w = Ua(4);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue


R(1) = q0_A * q0_A - (1. - q_A' * q_A);
R(5) = q0_B * q0_B - (1. - q_B' * q_B);

Ra(1:3) = p - (q0_A * q_B - q0_B * q_A - skew(q_A) * q_B);
Ra(4) = w * w - (1. - 0.25 * p' * p);

sigma = sys.parameters.K * p;

loc_gp{1} = [0.];
log_gp{2} = [];

weight_gp{1} = [2.];
weight_gp{2} = [];

n_gp = 1;
inda = 5;
for k = 1:n_gp
    
    alpha = 0.5 *(loc_gp{n_gp}(k) + 1.);
    
    zk = Ua(inda);
    bk = Ua(inda+1);
    fk = Ua(inda+2:inda+4);
    
    Ra(inda) = zk - (1. - 0.25 * alpha^2 * p' * p);
    Ra(inda+1) = zk * bk - p' * sigma;
    Ra(inda+2:inda+4) = zk * fk - (sigma + 0.25 * alpha^2 * bk * p);

    inda = inda + 5;
    
    wfk = w * fk;
    ck = 0.5 * skew(p) * fk;
    R(2:4) = R(2:4) + 0.5 * weight_gp{n_gp}(k) * ( -wfk - ck );
    R(6:8) = R(6:8) + 0.5 * weight_gp{n_gp}(k) * ( wfk - ck );
end

R(8) = R(8) + Moment;

R(2:4) = R(2:4) + Lambda_BC;
R(9:11) = q_A;

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
