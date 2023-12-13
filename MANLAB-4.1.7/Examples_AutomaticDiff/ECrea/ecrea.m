% from
%     Doedel, Keller, Kernevez (1991). 
%     Numerical analysis and control of bifurcation. 
%     Bifurcation and Chaos 1, 493-520.
%
%     Governing equations
%
%     2u1 - u2 + 100u1/(1+u1+u1*u1) = lambda 
%     2u2 - u1 + 100u2/(1+u2+u2*u2) = lambda+mu
%
%     Unknowns U = [ u1 u2 lambda ]  mu : fixed parameter

addpath('../../SRC')

%% Parameters of the system
neq = 2;
parameters.mu=0.1;

%% Initialization of the system :
sys = Syst('neq',neq,'parameters',parameters,'point_display',@point_display,'global_display',@global_display,'R',@R);

%% Starting point :
u0 = [0;0];
lambda0 = 0;
U0 = [u0; lambda0];

%% Manlab start :
Manlab('sys'            ,sys , ...      
       'U0value'         ,U0 , ...
       'displayvariables',[3 1; 3 2]);     % MANLAB run


