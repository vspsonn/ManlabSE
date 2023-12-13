%
% Bratu problem
% 
%
% Bratu continuous equations:
%   find U(x)=(u(x),lambda)) such that
%   u(x)'' + lambda exp(u(x)) = 0 in (0,1), u(0)=u(1)=0.
%
% Bratu discrete equations: 
% Interval (0,1) is divided into N+1 element -> N+2 node at 0 x1 x2   ... xN  1
% the unknowns are u1=u(x1),u2=u(x2), ... uN=u(xN)
% Laplacian is discretized using a central finite difference scheme
% 
%     -2u1 + u2                        + lambda exp(u1) =0
%       u1 -2u2 + u3                   + lambda exp(u2) =0
%                  ...
%                      u{N-2} -2u{N-1} + uN    + lambda exp(u{N-1}) =0
%                               u{N-1} -2uN    + lambda exp(uN)     =0
%
% or in vector form  A.u + lambda exp(u)=0
%
% The unknown vector is U=(u1,...uN,lambda)=(u,lambda)
%
%  Auxiliary variables : vi=exp(ui)   i=1..N
%
%  quadratic expression
%
%   A.u + lambda v =0
%   dv - Vdu =0 
%


global U Section Diagram    % Global variables to export from the diagram. 

addpath('../../SRC')

%% Parameters of the system
N=30;
neq = N;
neq_aux = N;

writing = 'standard';% Writing of the system {'standard' or 'vectorial'}
                     % The initialization is faster with vectorial form.

% Specific parameters
A=ones(N,1)*(N+1)^2; % Diagonal vector (with terms 1/dx^2)
K=spdiags([A -2*A A],-1:1,N,N); % Laplacien 1D DFC2 : centered scheme with 3 points (1 -2 1)/dx^2

parameters.N=N;
parameters.K=K;

%% initialization of the system
sys=SystAQ(neq,neq_aux,@equations_vec,@point_display,@global_display,parameters,writing);

%% starting point
lambda = 0; 
u0=zeros(N,1);

U0 = [u0 ; lambda ; exp(u0)];

%% Launch Manlab
Manlab('sys'            ,sys , ...      
         'U0value'         ,U0, ...
         'displayvariables',[neq+1 N/2]);     % MANLAB run

