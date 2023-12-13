% Pendulum : x'' + sin(x) = 0
% With unfolding term :
% x'' + lambda x' + sin(x) = 0
%
% Quadratic recast :
% x'' + lambda x' + s = 0
% s - sin(x)          = 0       | s' - x' * c = 0
% c - cos(x)          = 0       | c' + x' * s = 0

addpath(genpath('../../SRC'))

%% Parameters of the system
nz= 1;     % number of main equations of the DAE system
nz_aux =2;  % number of auxiliary equations of the DAE system
H = 30;     % number of harmonics
parameters = struct();

%% initialization of the system
type = 'conservative';        % type of system (can be 'forced' or 'autonomous' or 'conservative')
writing = 'standard';   % way the equations have been written (can be 'standard' or 'vectorial')

sys=SystHBQ(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,type,writing);

%% starting point
nz_tot=nz+nz_aux;
Z0  =ones(2*H+1,nz_tot)*1e-10;
lambda=0;omega=1;
Z0(2   ,1)=0.02 ;
% Auxiliary variables
Z0(:,2) = Z0(:,1);  % sin(x) ~= x for x<<1 
Z0(1,3) = 1;        % cos(x) ~= 1 for x<<1

%%% Initialization of the column vector of unknowns
U0 = sys.init_U0(Z0,omega,lambda);

% Display the first cosine and the first sine of the first variable with
% respect to the angular frequency omega.
dispvars = [sys.getcoord('omega') sys.getcoord('cos',1,1) ; ...
    sys.getcoord('omega') sys.getcoord('sin',1,1)];

%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,sys.init_U0(Z0,omega,lambda) , ...
    'ANMthreshold'    ,1e-16 , ...
    'NRthreshold'     ,1e-13, ...
    'NRmethod'        ,2, ...         % correction with variable parameter
    'displayvariables',dispvars);     % MANLAB run

