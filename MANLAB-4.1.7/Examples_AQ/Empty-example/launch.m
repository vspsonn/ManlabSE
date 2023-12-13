% Empty example.
% Write here the description of your system, the auxiliary variables used,
% and its final quadratic recast

global U Section Diagram   % Global variables to export point from the diagram. 

% Path of the SRC file.
addpath('../../SRC')

%% parameters of the system
neq = 2;     % Number of main equations
neq_aux = 6; % Number of auxiliary variables used

writing = 'standard';% Writing of the system {'standard' or 'vectorial'}

% Parameters specific to your system
parameters.param1 = 0.5;
parameters.param2 = 1;

%% initialization of the system - Should not be changed
sys=SystAQ(neq,neq_aux,@equations,@point_display,@global_display,parameters,writing);

%% starting point
%%% Main variables u and continuation parameter lambda
u = 0;
lambda = 5;

%%% Auxiliary variables v :
v = 3;

%%% Vector containing all the unknowns
U0 = [u lambda v]'; 

%% Launch Manlab
Manlab('sys'               ,sys, ...        % description of your system
         'U0value'         ,U0, ...         % starting point
         'displayvariables',[neq+1 1]);     % MANLAB run

% With 'displayvariables',[neq+1 1] ; the projected bifurcation diagram
% (figure 2) will represents the first variable with respect to the
% continuation parameter lambda.

%% Optional arguments of Manlab launching function with their default values :
%          'order'           ,20, ...     % order of the series
%          'ANMthreshold'    ,1e-6, ...   % threshold for the domain of validity of the series
%          'Amax_max'        ,1e6, ...    % maximum value of the domain of validity of the series
%          'NRthreshold'     ,2e-5, ...   % threshold for Newton-Raphson (NR) corrections
%          'NRitemax'        ,10, ...     % Maximum number of iteration of NR algorithm
%          'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
%          'NRmethod'        ,0, ...      % NR corrections on/[off]
%          'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
%          'PointDisplay'    ,1, ...      % Point display [on]/off
%          'GlobalDisplay'   ,1, ...      % Global display [on]/off
%          'StabilityCheck'  ,0, ...      % Stability computation on/[off]
%          'StabTol'         ,1e-6, ...   % Stability tolerance
%
%               The stability requires to write the system of equations in
%               the explicit ODE form X' = f(X). 
%               The residue function is then R(X) = f(X) = 0.
%               In all other cases, it gives wrong results.
