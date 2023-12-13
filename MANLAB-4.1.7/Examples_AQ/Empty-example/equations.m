function [Rf,dRf] = equations(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]   

u      = Uf(1:sys.neq);      % Main variables
lambda = Uf(sys.neq+1);      % Continuation parameter
Ua     = Uf(sys.neq+2:end);  % Auxiliary variables

du     = dUf(1:sys.neq);     % differential of Main variables
dlambda= dUf(sys.neq+1);     % differential of the continuation parameter
dUa    = dUf(sys.neq+2:end); % differential of Auxiliary variables

%% Parameters of the system

param1 = sys.parameters.param1;
param2 = sys.parameters.param2;

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue

% Equations of the main system          | differentiation of the non-quadratic equations
R(1) =  ;                                dR(1) = 
R(2) =  ;
...
    
% Definition of the auxiliary variables | differentiation of the non-quadratic equations
Ra(1) =  ;                             dRa(1) =
Ra(2) =  ;
...
    

% Concatenation of the two residues
Rf =[R  ; Ra ];
dRf=[dR ; dRa];

end
