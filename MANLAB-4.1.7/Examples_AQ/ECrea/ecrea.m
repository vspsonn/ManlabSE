%  ECrea example
%
%  Equations for ECREA 
%
%  r1 := 2 u1 - U2 + 100 u1/(1 + u1 + u1^2) = lambda
%  r2 := 2 u2 - u1 + 100 U2/(1 + u2 + u2^2) = lambda + mu
%
%  where u1,u2 are the unknowns and lambda the control parameter.
%  mu is a constant
%
%  By using the following auxiliary variables 
%
%  v1 = 1 + u1 + u1*u1
%  v2 = 1 + u2 + u2*u2
%  v3 = 1/v1
%  v4 = 1/v2 
%
%  the equation are written quadratic as
%
%   r1 = 2 u1 -  u2  + 100 u1 v3 -lambda       = 0
%   r2 = 2 u2 -  u1  + 100 u2 v4 -lambda - mu  = 0
%
%  The equation defining the auxiliary variables are also written quadratic as
%
%  raux1 = v1 - 1 - u1 - u1*u1  = 0
%  raux2 = v2 - 1 - u2 - u2*u2  = 0
%  raux3 = 1 - v3*v1            = 0
%  raux4 = 1 - v4*v2            = 0
%
%  Both systems are defined in "equation.m"
%
%  U   =[ u1 u2 lambda ]
%  Ua  =[ v1  v2  v3  v4 ]
%  Uf  =[ U Ua]

global U Section Diagram    % Global variables to export from the diagram. 

addpath('../../SRC')

%% parameters of the system
neq = 2;
neq_aux = 4;

% specific parameters
parameters.mu = 0.1; % There is a branching point if mu = 0

%% initialization of the system
sys=SystAQ(neq,neq_aux,@equations,@point_display,@global_display,parameters);
 
%% starting point
u1 = 0;
u2 = 0;
lambda = 0;

%%% Auxiliary variables
v1 = 1 + u1 + u1^2;
v2 = 1 + u2 + u2^2;
v3 = 1/v1;
v4 = 1/v2;

U0 = [u1 u2 lambda v1 v2 v3 v4]';

%% Launch Manlab
Manlab('sys'             ,sys, ...      
       'U0value'         ,U0, ...
       'displayvariables',[3 1; 3 2]);     % MANLAB run


