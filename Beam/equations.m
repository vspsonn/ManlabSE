function [Rf,dRf] = equations(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]

n_gp = sys.parameters.n_gp;

% n_nodes = 2;
% u0 = zeros(sys.neq+1,1);
% for i = 1:n_nodes
%    u0(1 + (i-1)*7) = 1; 
% end
% u0(12) = 1;
% 
% neq_aux = 14 + 21 * n_gp;
% ua0 = zeros(neq_aux,1);
% ua0(4) = 1;
% ua0(8) = 1;
% ua0(11) = 1;
% for i = 1:n_gp
%     ua0(15 + (i-1)*21) = 1;
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

Lambda_BC = U(15:20);
Moment = U(sys.neq+1);

p_theta = Ua(1:3);
w = Ua(4);

y = Ua(5:7);
d = Ua(8:10);
p_u = Ua(11:13);
rho = Ua(14);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue


R(1) = q0_A * q0_A - (1. - q_A' * q_A);
R(5) = q0_B * q0_B - (1. - q_B' * q_B);

Ra(1:3) = p_theta - (q0_A * q_B - q0_B * q_A - skew(q_A) * q_B);
Ra(4) = w * w - (1. - 0.25 * p_theta' * p_theta);

dx = x_B - x_A;
Ra(5:7) = y - skew(q_A) * dx;
Ra(8:10) = d - (dx - 2 * q0_A * y + 2 * (skew(q_A) * y));
Ra(11:13) = p_u - (w * d - 0.5 * (skew(p_theta) * d));
Ra(14) = rho - p_u' * p_theta;

f0 = zeros(6,1);
f0(1) = 1;
defo = [p_u - f0(1:3); p_theta - f0(4:6)];
sigma = sys.parameters.K * defo;

loc_gp{1} = [0.];
weight_gp{1} = [2.];

loc_gp{2} = [-1/sqrt(3), 1/sqrt(3)];
weight_gp{2} = [1., 1.];

inda = 15;
for k = 1:n_gp
    
    alpha = 0.5 *(loc_gp{n_gp}(k) + 1.);
    
    zk = Ua(inda);
    bk = Ua(inda+1);
    fk = Ua(inda+2:inda+4);
    
    Ra(inda) = zk - (1. - 0.25 * alpha^2 * p_theta' * p_theta);
    Ra(inda+1) = zk * bk - p_theta' * sigma(4:6);
    Ra(inda+2:inda+4) = zk * fk - (sigma(4:6) + 0.25 * alpha^2 * bk * p_theta);
    
    rhob_k = Ua(inda+5);
    pi_u = Ua(inda+6:inda+8);
    ksi_theta = Ua(inda+9);
    f_u = Ua(inda+10:inda+12);
    ksi_u = Ua(inda+13);
    ksib_u = Ua(inda+14);
    f_wu = Ua(inda+15:inda+17);
    fb_u = Ua(inda+18:inda+20);
    
    Ra(inda+5) = zk * rhob_k - rho;
    Ra(inda+6:inda+8) = pi_u - sys.parameters.K(1:3,1:3) * (p_u + 0.25 * alpha^2 * rhob_k * p_theta - w * f0(1:3));
    Ra(inda+9) = zk * ksi_theta - p_theta' * pi_u;
    Ra(inda+10:inda+12) = zk * f_u - (pi_u + 0.25 * alpha^2 * ksi_theta * p_theta);
    Ra(inda+13) = zk * ksi_u - p_u' * pi_u;
    Ra(inda+14) = ksib_u - (ksi_u + 0.75 * alpha^2 * rhob_k * ksi_theta);
    Ra(inda+15:inda+17) = zk * f_wu - 0.25 * alpha^2 * (p_u * ksi_theta + rhob_k * pi_u + ksib_u * p_theta);
    
    Ra(inda+18:inda+20) = w * fb_u - rho * f_u;
    
    inda = inda + 21;
    
    wfk = w * (fk + f_wu);
    wfu = w * f_u;
    ck = 0.5 * skew(p_theta) * (fk + f_wu);
    cu = 0.5 * skew(p_u) * f_u;
    ct = 0.5 * skew(p_theta) * f_u;
    R(2:4) = R(2:4) + 0.5 * weight_gp{n_gp}(k) * (0.25 * fb_u - cu - wfk - ck);
    R(5:7) = R(5:7) + 0.5 * weight_gp{n_gp}(k) * ( -wfu - ct );
    R(9:11) = R(9:11) + 0.5 * weight_gp{n_gp}(k) * (-0.25 * fb_u - cu + wfk - ck);
    R(12:14) = R(12:14) + 0.5 * weight_gp{n_gp}(k) * ( wfu - ct );
end

R(11) = R(11) + Moment;

R(2:7) = R(2:7) + Lambda_BC;
R(15:17) = q_A;
R(18:20) = x_A;

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
