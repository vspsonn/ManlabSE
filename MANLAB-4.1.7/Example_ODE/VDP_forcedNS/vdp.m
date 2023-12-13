% Launch Van der Pol problem
%
%    DAE for a forced van der pol model ( u' = u // v' - (mu1 - mu2 u - u^2) v + u = F cos(lambda t) )  (nz=2)
%
% The derivated terms should be with a -1 and the forced term is
% isolated on the right hand side, as follows :
%
%    v                    - u'  = 0
%  - u              + rv  - v'  = - F cos(lambda t)
%  mu1 - r - mu2*u  - u^2       = 0
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

global U Section    % Global variables to export point from the diagram.

addpath(genpath('../../SRC'))

%% Parameters of the system
nz= 2;     % number of main equations of the DAE system
nz_aux =1;  % number of auxiliary equations of the DAE system
H = 5*[1;1];     % number of harmonics

% specific parameters
parameters.mu1 = 0.1;
parameters.mu2 = 0.1;
parameters.F   = 1;

%% initialization of the system
sys=SystODE(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,'forced');

%% starting point

if numel(H) == 1
    %%% Periodic solution
    Z0  =randn(2*H+1,sys.nz_tot)*1e-7;
    omega = 0.3; lambda=omega;
    Z0(2,1) = .1;
    U0 = sys.init_U0(Z0,omega,lambda);
    
    % Display the first cosine and the first sine of the first variable with
    % respect to the bifurcation parameter lambda.
    dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,1) ; ...
        sys.getcoord('lambda') sys.getcoord('sin',1,1)];
    
    stabcheck = 1;
    
else
    %%% Quasi-periodic from a Neimark-Sacker bifurcation
    load('Section_NS_H10');
    U0 = sys.init_NS(Section);
    
    % Display the cosine (1,0) and the cosine (0,1) of the first variable with
    % respect to the bifurcation parameter lambda.
    dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,[1 0]) ; ...
        sys.getcoord('lambda') sys.getcoord('cos',2,[0 1])];
    
    stabcheck = 0;
    
end

%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'ANMthreshold'    ,1e-8 , ...
    'NRthreshold'     ,1e-5, ...
    'StabilityCheck'  ,stabcheck, ...
    'displayvariables',dispvars);     % MANLAB run


