% Launch Duffing problem
%
%    DAE for the forced duffing model ( u' = v // v' = - epsilon*v - u - u^3 = alpha*cos(lambda t)  (nz=2)
%    v                   - u'  = 0
%  - u - epsilon v - ru  - v'  = - alpha cos(lambda t)
%    r             - u^2       = 0
%
%    z(t) = [ u(t) ; v(t) ; r(t) ]  (nz_tot=3)  % unknown of the DAE
%
%    ^
%    u = [ u_0 ;  u_C1 ;  u_C2 ; ...  ; u_CH ; u_S1 ;  u_S2 ; ...  ; u_SH ]  Fourier coef (colum) vector for u
%
%          ^   ^     ^
%    Z= [ u , v , r ]  matrix of unknown Fourier coefficients
%
%    Z= reshape(Z, nz_tot*(2*H+1),1)  -> column vector of unknown Fourier coefficients

global U Section Diagram % Global variables to export point from the diagram. 

addpath(genpath('../../SRC'))

%% Parameters of the system
nz= 2;     % number of main equations of the DAE system
nz_aux =1;  % number of auxiliary equations of the DAE system
H = 25;     % number of harmonics

% specific parameters
parameters.alpha = 1;
parameters.epsilon = .05;

% Angular frequency. 
% Set to the prescribed value if this is not the continuation parameter.
parameters.angfreq = 'omega';

%% initialization of the system
type = 'forced';        % type of system (can be 'forced' or 'autonomous')
writing = 'standard';   % way the equations have been written (can be 'standard' or 'vectorial')

sys=SystHBQ(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,type,writing);

%% starting point
omega=0.4;lambda=omega;

Z0  =randn(2*H+1,sys.nz_tot)*1e-7;
Z0(2   ,1)=0.7;
Z0(:,2)=omega*sys.D(Z0(:,1));        % v = u' in the Fourier domain.
Z0(:,3) = sys.Prod(Z0(:,1),Z0(:,1)); % r = u*u in the Fourier domain.

U0 = sys.init_U0(Z0,omega,lambda);

% Display the first cosine and the first sine of the first variable with
% respect to the bifurcation parameter lambda.
dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,1) ; ...
    sys.getcoord('lambda') sys.getcoord('sin',1,1)];

%% Launch of Manlab with options
% Manlab('sys'       ,sys , ...
%     'U0value'         ,U0, ...
%     'ANMthreshold'    ,1e-10 , ...
%     'NRthreshold'     ,1e-6, ...
%     'StabilityCheck'  ,1, ...
%     'StabTol'         ,1e-3, ...
%     'displayvariables',dispvars);     % MANLAB run
% 

 %% Manlab without interface
     
Diagram = Manlab_script('nb_step', 70 , ... % 70 continuation steps
    'sys'             ,sys , ...
    'U0value'         ,U0, ...
    'ANMthreshold'    ,1e-10 , ...
    'NRthreshold'     ,1e-6, ...
    'StabilityCheck'  ,1, ...
    'StabTol'         ,1e-3, ...
    'displayvariables',dispvars);     % MANLAB run

 
%%% Post processing
global_display(sys,Diagram);

