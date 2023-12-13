% Launch exp_tanh example. 
% Guillot Cochelin Vergez 2018 "A generic and efficient Taylor series based
% continuation method using a quadratic recast of smooth nonlinear systems"
% https://hal.archives-ouvertes.fr/hal-01827832
%
% Equations for exp-tanh example :
%  u1 + lambda*exp(u2)/(1+u1)   = 0
%  u2 + u1*tanh(-5u1/(1+u1*u2)) = 0
%
%   Can be written quadratically with the help of 6 auxiliary variables.
%
%   v1 = 1 + u1*u2;
%   v2 = -5*u1/v1;
%   v3 = tanh(v2);
%   v4 = 1-v3^2;
%   v5 = exp(u2);
%   v6 = v5/(1+u1);
%
%
%   Quadratic recast :
%
%   R(1) = 0   +   u1   +  lambda*v6
%   R(2) = 0   +   u2   +  u1*v3
%  Ra(1) =-1   +   v1   -  u1*u2
%  Ra(2) = 0   +  5u1   +  v1*v2
%  Ra(3) = 0   +   v3   +              -   tanh(v2)
%  Ra(4) = -1   +   v4   +  v3*v3
%  Ra(5) = 0   +   v5   +              -   exp (u2)
%  Ra(6) = 0   +  v6-v5 +  u1*v6
%
%   Rf = [R ; Ra]
%   Uf = [U lambda Ua] = [ u1 u2 lambda v1 v2 v3 v4 v5 v6]


global U Section Diagram    % Global variables to export from the diagram. 

addpath('../../SRC')

%% parameters of the system
neq = 2;
neq_aux = 6;

parameters = struct(); % no specific parameters.

%% initialization of the system
sys=SystAQ(neq,neq_aux,@equations,@point_display,@global_display,parameters);

%% starting point
u1 = 0;
u2 = 0;
lambda = 0;

%%% Auxiliary variables :
v1 = 1 + u1*u2;
v2 = -5*u1/v1;
v3 = tanh(v2);
v4 = 1-v3^2;
v5 = exp(u2);
v6 = v5/(1+u1);

U0 = [u1 u2 lambda v1 v2 v3 v4 v5 v6]'; 

%% Launch Manlab
Manlab('sys'            ,sys , ...      
         'U0value'         , U0, ...
         'displayvariables',[3 1; 3 2]);     % MANLAB run

