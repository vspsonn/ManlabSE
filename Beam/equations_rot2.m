function [Rf,dRf] = equations_rot2(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

loc_gp{1} = [0.];
weight_gp{1} = [2.];

loc_gp{2} = [-1/sqrt(3), 1/sqrt(3)];
weight_gp{2} = [1., 1.];

loc_gp{3} = [0.0, -sqrt(3/5), sqrt(3/5)];
weight_gp{3} = [8.0/9.0, 5.0/9.0, 5.0/9.0];


%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]

n_gp = sys.parameters.n_gp;

% n_nodes = 2;
% u0 = zeros(sys.neq+1,1);
% for i = 1:n_nodes
%     u0(1 + (i-1)*4) = 1;
% end
% ua0 = zeros(sys.neq_aux,1);
% inda = 1;
% for k = 1:n_gp
%     ua0(inda) = 1;
%     ua0(inda+4) = 1;
%     ua0(inda+8) = 1;
%     inda = inda + (4 * n_nodes + 4);
% end
% Uf = [u0; ua0];

U = Uf(1:sys.neq+1);      % Main variables
Ua = Uf(sys.neq+2:end);  % Auxiliary variables

q0_A = U(1);
q_A = U(2:4);
q0_B = U(5);
q_B = U(6:8);

Lambda_BC = U(9:11);
Moment = U(sys.neq+1);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue

R(1) = q0_A * q0_A - (1. - q_A' * q_A);
R(5) = q0_B * q0_B - (1. - q_B' * q_B);

inda = 1;
Q = zeros(3,6);
kappa = zeros(3,1);
for k = 1:n_gp
    
    J0 = 0.5;
    w_gp = weight_gp{n_gp}(k) * J0;
    
    eta = loc_gp{n_gp}(k);
    
    N_A = 0.5*(1-eta);
    N_B = 0.5*(1+eta);
    dN_A = -0.5;
    dN_B = 0.5;
    
    qb = N_A * [q0_A; q_A] + N_B * [q0_B; q_B];
    qb0 = qb(1);
    qb_ = qb(2:4);
    
    zk = Ua(inda);
    Ra(inda) = zk * zk - qb' * qb;
    inda = inda + 1;
    
    pk_A = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk * pk_A + 2.0 * (q0_A * qb_ - qb0 * q_A - skew(q_A) * qb_);
    inda = inda + 3;
   
    wk_A = Ua(inda);
    Ra(inda) = wk_A * wk_A - (1.0 - 0.25 * pk_A' * pk_A);
    inda = inda + 1;
    
    pk_B = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk * pk_B + 2.0 * (q0_B * qb_ - qb0 * q_B - skew(q_B) * qb_);
    inda = inda + 3;
    
    wk_B = Ua(inda);
    Ra(inda) = wk_B * wk_B - (1.0 - 0.25 * pk_B' * pk_B);
    inda = inda + 1;
    
    wb = N_A * wk_A + N_B * wk_B;
    
    yk = Ua(inda:inda+2);
    Ra(inda:inda+2) = wb * yk - (dN_A/J0 * pk_A + dN_B/J0 * pk_B);
    inda = inda + 3;
    
%     kappa = kappa + w_gp * yk;
%     Q = Q + w_gp * [dN_A/J0*eye(3)+N_A*skew(yk) dN_B/ J0*eye(3)+N_B*skew(yk)];
    
        R([2:4 6:8]) = R([2:4 6:8]) ...
            + [dN_A/J0*eye(3)+N_A*skew(yk) dN_B/J0*eye(3)+N_B*skew(yk)]' * (w_gp * sys.parameters.K * yk);
end

% R([2:4 6:8]) = Q' * (sys.parameters.K * kappa);
R(8) = R(8) + Moment;

R(2:4) = R(2:4) + Lambda_BC;
R(9:11) = q_A;

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
