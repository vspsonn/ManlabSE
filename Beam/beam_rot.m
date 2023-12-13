% Continuation of a simple exemple: the union of two parabolas:
% 
% first parabola: x=y^2
% second parabola: y=(x-2)^2+a where a is a parameter
% 
% Algebraic equation for the system:
% 
% r(x,y) = (x-y^2)(y-(x-2)^2+a) = 0
%
% Auxilary variables to obtain a quadratic system:
%
% u1 = x - y^2;
% u2 = y - (x-2)^2 + a
%
% Primary equations:
%
% r = u1 * u2 = 0
%
% Auxilary equations:
%
% raux1 = u1 -x + y^2
% raux2 = u2 -y + (x-2)^2 + a
%
% U   = [x y]
% Uaux = [u1 u2]
% Uf  = [U Uaux]


global U Section Diagram   % Global variables to export point from the diagram. 

% Path of the SRC file.
addpath('../../SRC')

%% parameters of the system

n_nodes = 2;
n_elem = n_nodes-1;
n_gp = 1;
n_BC = 3;

neq = n_BC + n_nodes * 4;     % Number of main equations
neq_aux = n_elem * 4 + (n_elem * n_gp * 5); % Number of auxiliary variables used

writing = 'standard'; % Writing of the system {'standard' or 'vectorial'}

% Parameters specific to your system
K=eye(3);
parameters.K = K;

%% initialization of the system - Should not be changed
sys=SystAQ(neq,neq_aux,@equations_rot,@point_display,@global_display,parameters,writing);

%% starting point
u0 = zeros(neq+1,1);
for i = 1:n_nodes
   u0(1 + (i-1)*4) = 1; 
end
ua0 = zeros(neq_aux,1);
inda = 5;
ua0(4) = 1;
for k = 1:n_gp
    ua0(inda) = 1;
    inda = inda + 5;
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
