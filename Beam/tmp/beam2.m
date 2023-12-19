global U Section Diagram   % Global variables to export point from the diagram.

% Path of the SRC file.
addpath('../../SRC')

loc_gp{1} = [0.];
weight_gp{1} = [2.];

loc_gp{2} = [-1/sqrt(3), 1/sqrt(3)];
weight_gp{2} = [1., 1.];

loc_gp{3} = [0.0, -sqrt(3/5), sqrt(3/5)];
weight_gp{3} = [8.0/9.0, 5.0/9.0, 5.0/9.0];

%% parameters of the system

n_nodes = 2;
n_elem = n_nodes-1;
n_gp = 2;
n_BC = 6;

neq = n_BC + n_nodes * 7 + 3 * n_gp;     % Number of main equations
neq_aux = n_elem * n_gp * (4 * n_nodes + 4 + 11 * n_nodes + 3) + 0; % Number of auxiliary variables used

writing = 'standard'; % Writing of the system {'standard' or 'vectorial'}

% Parameters specific to your system
K=eye(6);
parameters.K = K;
parameters.n_gp = n_gp;
parameters.n_nodes = n_nodes;

%% initialization of the system - Should not be changed
sys=SystAQ(neq,neq_aux,@equations2,@point_display,@global_display,parameters,writing);

%% starting point
u0 = zeros(neq+1,1);
for i = 1:n_nodes
    u0(1 + (i-1)*7) = 1;
end
u0(12) = 1;

q0_A = u0(1);
q_A = u0(2:4);
x_A = u0(5:7);
q0_B = u0(8);
q_B = u0(9:11);
x_B = u0(12:14);

ua0 = zeros(neq_aux,1);
inda = 1;
ind = 15;
for k = 1:n_gp
    J0 = 0.5;
    eta = loc_gp{n_gp}(k);
    
    N_A = 0.5*(1-eta);
    N_B = 0.5*(1+eta);
    dN_A = -0.5;
    dN_B = 0.5;
    
    qb = N_A * [q0_A; q_A] + N_B * [q0_B; q_B];
    qb0 = qb(1);
    qb_ = qb(2:4);
    
    zk = sqrt(qb' * qb);
    ua0(inda) = zk;
    inda = inda + 1;
    
    pk_A = - 2.0 / zk * (q0_A * qb_ - qb0 * q_A - skew(q_A) * qb_) ;
    ua0(inda:inda+2) = pk_A;
    inda = inda + 3;
    
    pk_B = -2.0 / zk * (q0_B * qb_ - qb0 * q_B - skew(q_B) * qb_);
    ua0(inda:inda+2) = pk_B;
    inda = inda + 3;
    
    wk_A = sqrt(1.0 - 0.25 * (pk_A' * pk_A));
    ua0(inda) = wk_A;
    inda = inda + 1;
    
    wk_B = sqrt(1.0 - 0.25 * (pk_B' * pk_B));
    ua0(inda) = wk_B;
    inda = inda + 1;
    
    wb = N_A * wk_A + N_B * wk_B;
    
    yk = (dN_A/J0 * pk_A + dN_B/J0 * pk_B) / wb;
    ua0(inda:inda+2) = yk;
    inda = inda + 3;
    
    xk = [0.5 *(loc_gp{n_gp}(k) + 1.); 0; 0];
    u0(ind:ind+2) = xk;
    ind = ind + 3;
    
    uk_A = skew(qb_) * (x_A - xk) / zk;
    ua0(inda:inda+2) = uk_A;
    inda = inda + 3;
    
    uk_B = skew(qb_) * (x_B - xk) / zk;
    ua0(inda:inda+2) = uk_B;
    inda = inda + 3;
    
    dk_A = (x_A - xk) - 2.0 / zk * (qb0 * uk_A - skew(qb_) * uk_A);
    ua0(inda:inda+2) = dk_A;
    inda = inda + 3;
    
    dk_B = (x_B - xk) - 2.0 / zk * (qb0 * uk_B - skew(qb_) * uk_B);
    ua0(inda:inda+2) = dk_B;
    inda = inda + 3;
    
    puk_A = wk_A * dk_A - 0.5 * skew(pk_A) * dk_A;
    ua0(inda:inda+2) = puk_A;
    inda = inda + 3;
    
    puk_B = wk_B * dk_B - 0.5 * skew(pk_B) * dk_B;
    ua0(inda:inda+2) = puk_B;
    inda = inda + 3;
    
    bk_A = - 0.25 * puk_A' * pk_A / wk_A;
    ua0(inda) = bk_A;
    inda = inda + 1;
    
    bk_B = - 0.25 * puk_B' * pk_B / wk_B;
    ua0(inda) = bk_B;
    inda = inda + 1;
    
    rhob = N_A * bk_A + N_B * bk_B;
    
    yuk = ((dN_A/J0 * puk_A + dN_B/J0 * puk_B) - rhob * yk) / wb;
    ua0(inda:inda+2) = yuk;
    inda = inda + 3;
end
U0 = [u0; ua0];

%% Launch Manlab
Manlab('sys'               ,sys, ...        % description of your system
    'U0value'         ,U0, ...         % starting point
    'displayvariables',[neq+1 8], ...    % MANLAB run
    'order'           ,20, ...     % order of the series
    'ANMthreshold'    ,1e-8, ...   % threshold for the domain of validity of the series
    'Amax_max'        ,1e6, ...    % maximum value of the domain of validity of the series
    'NRthreshold'     ,1e-7, ...   % threshold for Newton-Raphson (NR) corrections
    'NRitemax'        ,10, ...     % Maximum number of iteration of NR algorithm
    'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
    'NRmethod'        ,0, ...      % NR corrections on/[off]
    'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
    'PointDisplay'    ,0, ...      % Point display [on]/off
    'GlobalDisplay'   ,0, ...      % Global display [on]/off
    'StabilityCheck'  ,0, ...      % Stability computation on/[off]
    'StabTol'         ,1e-6);   % Stability tolerance

% With 'displayvariables',[neq+1 1] ; the projected bifurcation diagram
% (figure 2) will represents the first variable with respect to the
% continuation parameter lambda.
%
%               The stability requires to write the system of equations in
%               the explicit ODE form X' = f(X).
%               The residue function is then R(X) = f(X) = 0.
%               In all other cases, it gives wrong results.
