%
% Logistic map
% x_{n+1} = f(x_n) := mu x_n (1-x_n)
%
% Continue the N-th first iterates (x_1,x_2,x_3,...,x_N) of the sequence
% with respect to mu.
%
% It undergoes a large number of flip bifurcations when mu increases toward
% the critical value of 3.57.. for which the system becomes chaotic.


global U Section Diagram   % Global variables to export point from the diagram.

addpath('../../SRC')

%% Parameters of the system
N=5000; % number of points of the sequence
neq = N;
neq_aux = 1;

% Specific parameters
parameters.N=N;
parameters.x0=0.5;

%% initialization of the system
eq = @equations_vec;
writing = 'vectorial';
%%% Generation of the operators C,L,Q of the system by hand. This allows a
%%% much mor efficient initialization of the system.
operators = get_operators(parameters);

sys=SystAQ(neq,neq_aux,eq,@point_display,@global_display,parameters,writing,operators);

%% starting point
mu = 0.01; % in [0,4]
x = zeros(N,1);
x(1)=parameters.x0;
for i=2:N
    x(i) = mu*x(i-1)*(1-x(i-1));
end

U0 = [x;mu;1/mu];

dispvars = [ones(8,1)*(neq+1) , (N-7:N)'];

%% Launch Manlab
Manlab('sys'             ,sys , ...
    'NRthreshold'     ,1e-4, ...
    'ANMthreshold'    ,1e-6, ...
    'U0value'         ,U0, ...
    'GlobalDisplay'   ,0, ...
    'displayvariables',dispvars);     % MANLAB run