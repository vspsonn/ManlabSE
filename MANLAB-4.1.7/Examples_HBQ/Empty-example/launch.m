% Empty example.
% Write here the description of your system, the auxiliary variables used,
% and its final quadratic recast

global U Section Diagram   % Global variables to export point from the diagram.

% Path of the SRC file.
addpath(genpath('../../SRC'));

%% Parameters of the system
nz= 1;               % number of main equations of the system of the differential-algebraic system (DAE)
nz_aux =2;           % number of auxiliary equations of the system of the DAE

H = 20;              % number of harmonics used to compute the solution-branch

%%% Parameters specific to your system
parameters.param1 = 0.5;
parameters.param2 = 1;

parameters.angfreq=1; % if the system is forced at a fixed angular frequency, parameters.angfreq
% contains its value. Otherwise, parameters.angfreq = 'omega'.

%% initialization of the system
type = 'forced';        % type of system (can be 'forced' or 'autonomous')
writing = 'standard';   % way the equations have been written (can be 'standard' or 'vectorial')

sys=SystHBQ(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,type,writing);
%% starting point
lambda0=0;   % Continuation parameter initial value
omega0=1;    % Pulsation initial value

Z0 =zeros(2*H+1,sys.nz_tot);
% The matrix Z0 contains in column the Fourier development of all the
% variables of your system :
% Z0 = [ u , v ]
% where u = [ u_0 ; u_c1 ; u_c2 ; ... ; u_cH ; u_s1 ; u_s2 ; ... ; u_sH ];
% u_0 is the constant Fourier coefficient, u_c1 the first cosine Fourier
% coefficient, u_s1 the first sine Fourier coefficient, etc...
Z0(1,1)     = u0;
Z0(2,1)     = uc1;
Z0(1+j,1)   = ucj; % with 1 <= j <= H
Z0(2+H,1)   = us1;
Z0(1+H+j,1) = usj; % with 1 <= j <= H

%%% Vector U0 containing the all the unknowns.
U0 = sys.init_U0(Z0,omega0,lambda0);

%%% Variable displayed in the projected bifurcation diagram.
% To plot the coefficient of cos(h omega t) of variable number i with
% respect to lambda you should write as follow :
dispvars = [sys.getcoord('lambda') sys.getcoord('cos',i,h)];

%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'displayvariables',dispvars);     % MANLAB run


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
