%  Equations for VDP (x'' + alpha x' - lambda(1-x^2)x' + x = 0) :
%
% Order 1 recast :
%  x' = y
%  y' = -alpha y + lambda (1 - x^2)y - x
%
% Two auxiliary variables are added :
%   v = lambda y
%   r = 1-x^2
%
% Quadratic recast :
%  x' = y
%  y' = -alpha y + r*v - x
%
% This system is known to encounter a Hopf bifurcation for lambda = alpha.

global U Section Diagram    % Global variables to export point from the diagram.

addpath('../../SRC')

%% parameters of the system
nz = 2;
nz_aux = 2;
H = 5; % or H=0 for equilibrium continuation

% specific parameters
parameters.alpha = 0.5;

%% initialization of the system
sys=SystODE(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,'autonomous');

%% starting point

if H == 0
    x = 0;
    y = 0;
    lambda = 0;
    
    %%% Auxiliary variables
    v = lambda*y;
    r = 1-x^2;
    
    U0 = [x y lambda v r]';
    
    dispvars = [3 1; 3 2];
    
else
    
    load('Section_Hopf');
    U0 = sys.init_Hopf(Section);
    
    dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,1) ; ...
        sys.getcoord('lambda') sys.getcoord('omega')];
end

%% Launch Manlab
Manlab('sys'             ,sys, ...
    'U0value'         ,U0, ...
    'NRthreshold'     ,1e-10, ...
    'ANMthreshold'    ,1e-14, ...
    'Amax_max'        ,.7, ...     % Maximum value of the domain of validity of the series.
    'StabilityCheck'  ,1, ...
    'displayvariables',dispvars);     % MANLAB run


