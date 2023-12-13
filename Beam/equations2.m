function [Rf,dRf] = equations2(sys,Uf,dUf)
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
n_nodes = sys.parameters.n_nodes;

% u0 = zeros(sys.neq+1,1);
% for i = 1:n_nodes
%     u0(1 + (i-1)*7) = 1;
% end
% u0(12) = 1;
% 
% q0_A = u0(1);
% q_A = u0(2:4);
% x_A = u0(5:7);
% q0_B = u0(8);
% q_B = u0(9:11);
% x_B = u0(12:14);
% 
% ua0 = zeros(sys.neq_aux,1);
% inda = 1;
% ind = 15;
% for k = 1:n_gp
%     J0 = 0.5;
%     eta = loc_gp{n_gp}(k);
%     
%     N_A = 0.5*(1-eta);
%     N_B = 0.5*(1+eta);
%     dN_A = -0.5;
%     dN_B = 0.5;
%     
%     qb = N_A * [q0_A; q_A] + N_B * [q0_B; q_B];
%     qb0 = qb(1);
%     qb_ = qb(2:4);
%     
%     zk = sqrt(qb' * qb);
%     ua0(inda) = zk;
%     inda = inda + 1;
%     
%     pk_A = - 2.0 / zk * (q0_A * qb_ - qb0 * q_A - skew(q_A) * qb_) ;
%     ua0(inda:inda+2) = pk_A;
%     inda = inda + 3;
%     
%     pk_B = -2.0 / zk * (q0_B * qb_ - qb0 * q_B - skew(q_B) * qb_);
%     ua0(inda:inda+2) = pk_B;
%     inda = inda + 3;
%     
%     wk_A = sqrt(1.0 - 0.25 * (pk_A' * pk_A));
%     ua0(inda) = wk_A;
%     inda = inda + 1;
%     
%     wk_B = sqrt(1.0 - 0.25 * (pk_B' * pk_B));
%     ua0(inda) = wk_B;
%     inda = inda + 1;
%     
%     wb = N_A * wk_A + N_B * wk_B;
%     
%     yk = (dN_A/J0 * pk_A + dN_B/J0 * pk_B) / wb;
%     ua0(inda:inda+2) = yk;
%     inda = inda + 3;
%     
%     xk = [0.5 *(loc_gp{n_gp}(k) + 1.); 0; 0];
%     u0(ind:ind+2) = xk;
%     ind = ind + 3;
%     
%     uk_A = skew(qb_) * (x_A-xk) / zk;
%     ua0(inda:inda+2) = uk_A;
%     inda = inda + 3;
%     
%     uk_B = skew(qb_) * (x_B-xk) / zk;
%     ua0(inda:inda+2) = uk_B;
%     inda = inda + 3;
%     
%     dk_A = (x_A-xk) - 2.0 / zk * (qb0*uk_A - skew(qb_) * uk_A);
%     ua0(inda:inda+2) = dk_A;
%     inda = inda + 3;
%     
%     dk_B = (x_B-xk) - 2.0 / zk * (qb0*uk_B - skew(qb_) * uk_B);
%     ua0(inda:inda+2) = dk_B;
%     inda = inda + 3;
%     
%     puk_A = wk_A * dk_A - 0.5 * skew(pk_A) * dk_A;
%     ua0(inda:inda+2) = puk_A;
%     inda = inda + 3;
%     
%     puk_B = wk_B * dk_B - 0.5 * skew(pk_B) * dk_B;
%     ua0(inda:inda+2) = puk_B;
%     inda = inda + 3;
%     
%     bk_A = - 0.25 * puk_A' * pk_A / wk_A;
%     ua0(inda) = bk_A;
%     inda = inda + 1;
%     
%     bk_B = - 0.25 * puk_B' * pk_B / wk_B;
%     ua0(inda) = bk_B;
%     inda = inda + 1;
%     
%     rhob = N_A * bk_A + N_B * bk_B;
%     
%     yuk = -(rhob * yk - (dN_A/J0 * puk_A + dN_B/J0 * puk_B)) / wb;
%     ua0(inda:inda+2) = yuk;
%     inda = inda + 3;
% end
% Uf = [u0; ua0];

U = Uf(1:sys.neq+1);      % Main variables
Ua = Uf(sys.neq+2:end);  % Auxiliary variables

q0_A = U(1);
q_A = U(2:4);
x_A = U(5:7);
q0_B = U(8);
q_B = U(9:11);
x_B = U(12:14);

Lambda_BC = U(sys.neq-5:sys.neq);
Moment = U(sys.neq+1);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue

R(1) = q0_A * q0_A - (1. - q_A' * q_A);
R(5) = q0_B * q0_B - (1. - q_B' * q_B);

inda = 1;
ind = 15;
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
    
    pk_B = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk * pk_B + 2.0 * (q0_B * qb_ - qb0 * q_B - skew(q_B) * qb_);
    inda = inda + 3;
    
    wk_A = Ua(inda);
    Ra(inda) = wk_A * wk_A - (1.0 - 0.25 * (pk_A' * pk_A));
    inda = inda + 1;
    
    wk_B = Ua(inda);
    Ra(inda) = wk_B * wk_B - (1.0 - 0.25 * (pk_B' * pk_B));
    inda = inda + 1;
    
    wb = N_A * wk_A + N_B * wk_B;
    
    yk = Ua(inda:inda+2);
    Ra(inda:inda+2) = wb * yk - (dN_A/J0 * pk_A + dN_B/J0 * pk_B);
    inda = inda + 3;
    
    xk = U(ind:ind+2);
    
    uk_A = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk*uk_A - skew(qb_) * (x_A-xk);
    inda = inda + 3;
    
    uk_B = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk*uk_B - skew(qb_) * (x_B-xk);
    inda = inda + 3;
    
    dk_A = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk*dk_A - (zk*(x_A-xk) - 2.0 * (qb0 * uk_A - skew(qb_) * uk_A));
    inda = inda + 3;
    
    dk_B = Ua(inda:inda+2);
    Ra(inda:inda+2) = zk*dk_B - (zk*(x_B-xk) - 2.0 * (qb0 * uk_B - skew(qb_) * uk_B));
    inda = inda + 3;
    
    puk_A = Ua(inda:inda+2);
    Ra(inda:inda+2) = puk_A - (wk_A * dk_A - 0.5 * skew(pk_A) * dk_A);
    inda = inda + 3;
    
    puk_B = Ua(inda:inda+2);
    Ra(inda:inda+2) = puk_B - (wk_B * dk_B - 0.5 * skew(pk_B) * dk_B);
    inda = inda + 3;
    
    R(ind:ind+2) = N_A * puk_A + N_B * puk_B;
    ind = ind + 3;
    
    bk_A = Ua(inda);
    Ra(inda) = 4.0 * wk_A * bk_A + puk_A' * pk_A;
    inda = inda + 1;
    
    bk_B = Ua(inda);
    Ra(inda) = 4.0 * wk_B * bk_B + puk_B' * pk_B;
    inda = inda + 1;
    
    rhob = N_A * bk_A + N_B * bk_B;
    
    yuk = Ua(inda:inda+2);
    Ra(inda:inda+2) = wb * yuk + rhob * yk - (dN_A/J0 * puk_A + dN_B/J0 * puk_B);
    inda = inda + 3;
    
    defo = [yuk; yk] - [1;zeros(5,1)];
    stress =  sys.parameters.K * defo;
    
    yk_hat = [skew(yk) skew(yuk); zeros(3) skew(yk)];
    Q = [dN_A/J0*eye(6)+N_A*yk_hat dN_B/J0*eye(6)+N_B*yk_hat];
    
    R([2:7 9:14]) = R([2:7 9:14]) + Q' * (w_gp * stress);
end

R(11) = R(11) + Moment;

R(2:7) = R(2:7) + Lambda_BC;
R(ind:ind+5) = [q_A; x_A];

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
