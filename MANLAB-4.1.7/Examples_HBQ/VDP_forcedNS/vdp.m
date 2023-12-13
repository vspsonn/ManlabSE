% Launch Van der Pol problem, 
% see Guillot Vigué Vergez Cochelin 2017 "Continuation of quasi-periodic
% solutions with two-frequency harmonic balance method"
%
%    DAE for a forced van der pol model ( u' = u // v' + (mu1 + mu2 u - u^2) v + u = F cos(lambda t) )  (nz=2)
%
% The derivated terms should be with a -1 sign and the forced term is
% isolated on the right hand side, as follows :
%
%    v                    - u'  = 0
%  - u              - rv  - v'  = - F cos(lambda t)
%  mu1 - r - mu2*u  - u^2       = 0
%
%    z(t) = [ u(t) ; v(t) ; r(t) ]  (nz_tot=3)  % unknown of the DAE
%
%    U = [ u_0 ;  u_C1 ;  u_C2 ; ...  ; u_CH ; u_S1 ;  u_S2 ; ...  ; u_SH ]  Fourier coef (colum) vector for u
%
%    Z= [ U , V , R ]  matrix of unknown Fourier coefficients
%

global U Section Diagram    % Global variables to export from the diagram. 

addpath(genpath('../../SRC'))

%% Parameters of the system
nz= 2;     % number of main equations of the DAE system
nz_aux =1;  % number of auxiliary equations of the DAE system
H = 15;     % number of harmonics

% specific parameters
parameters.mu1 = 0.2;
parameters.mu2 = 0.1;
parameters.F   = 1;

% Angular frequency. 
parameters.angfreq = 'omega'; % if the system is forced at a fixed angular frequency, parameters.angfreq
                              % contains its value. Otherwise, parameters.angfreq = 'omega'.

%% initialization of the system
type = 'forced';        % type of system (can be 'forced' or 'autonomous')
writing = 'standard';   % way the equations have been written (can be 'standard' or 'vectorial')

sys=SystHBQ(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,type,writing);

%% starting point
omega = 0.6; lambda=omega;

Z0  =randn(2*H+1,sys.nz_tot)*1e-7;
Z0(2,1) = 0.5;
%%% Auxiliary variables :
Z0(:,2) = omega*sys.D(Z0(:,1));
Z0(1,3) = -parameters.mu1;
Z0(:,3) = Z0(:,3) + parameters.mu2*Z0(:,1) + sys.Prod(Z0(:,1),Z0(:,1));

%%% Initialization of the column vector of unknowns
U0 = sys.init_U0(Z0,omega,lambda);

% Display the first cosine and the first sine of the first variable with
% respect to the bifurcation parameter lambda.
dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,1) ; ...
    sys.getcoord('lambda') sys.getcoord('sin',1,1)];



interface = 1;

if interface == 1
%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'ANMthreshold'    ,1e-7 , ...
    'NRthreshold'     ,1e-5, ...
    'StabilityCheck'  ,1, ...
    'displayvariables',dispvars);     % MANLAB run


else
%% Manlab without interface

Diagram = Manlab_script('nb_step', 15 , ... % 15 continuation steps
    'sys'             ,sys , ...
    'U0value'         ,U0, ...
    'ANMthreshold'    ,1e-7 , ...
    'NRthreshold'     ,1e-5, ...
    'StabilityCheck'  ,0, ...               % Without stability
    'displayvariables',dispvars);     % MANLAB run

%%% Post - processing
% L2 norm of the first variable plotted versus lambda in figure 6.
plotdiagnormHBMbif(sys,Diagram,1,'lambda',6); 

% Concatenante all the solution-point in Diag.DiagUpp
Diag = calcdiagUpp(sys,Diagram);

% Values of lambda/omega along the diagram
lambda_diag = Diag.DiagUpp(sys.neq+1,:);
omega_diag = Diag.DiagUpp(sys.neq,:);

% Floquet exponents and multiplers computation :
% [Floquet_exponents,Floquet_multipliers] = calcdiagStab(sys,Diagram);

%%% Parallel post-processing for stability
[Floquet_exponents,Floquet_multipliers] = calcdiagStab_parallel(sys,Diagram);


figure(5)
subplot(2,2,[2 4]) 
plot(Floquet_multipliers.','o');hold on;
plot(cos(0:.01:2*pi),sin(0:.01:2*pi),'k'); % circle
subplot(2,2,1)
plot(lambda_diag,abs(Floquet_multipliers),'o')
subplot(2,2,3)
plot(lambda_diag,angle(Floquet_multipliers),'o')

end