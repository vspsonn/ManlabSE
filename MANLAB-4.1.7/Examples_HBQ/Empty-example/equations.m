function [Rf,dRa,Forcing] = equations(sys,t,zf,dzf,d2zf)
% Equations of the system of the form 
% Rf(zf) = C + L0(zf) + lambda L1(zf) + D0(dzf) + lambda D1(dzf) + DD(d2zf) + Q(zf,zf) + f(zf).

%%% parameters of the system
param1 = sys.parameters.param1;
param2 = sys.parameters.param2;

%% Variables : 

u      = zf(1:sys.nz);         % Main variables
Ua     = zf(sys.nz+1:end-1);   % Auxiliary variables
lambda = zf(end);              % Continuation parameter

du     = dzf(1:sys.nz);        % first order derivated of Main variables
dUa    = dzf(sys.nz+1:end);    % first order derivated of Auxiliary variables

d2u    = d2zf(1:sys.nz);       % second order derivated of Main variables
d2Ua   = d2zf(sys.nz+1:end);   % first order derivated of Auxiliary variables

%% Residues
R     = zeros(sys.nz,1);         % Main residue
Ra  = zeros(sys.nz_aux,1);     % Auxiliary residue

dRa = zeros(sys.nz_aux,1);     % Differential form of non-quadratic part of the auxiliary residue

% Equations of the main system
R(1) =  ;


% Definition of the auxiliary variables | differentiation of the non-quadratic equations
Ra(1) =  ;                             dRa(1) =

    
% Concatenation of the two residues
Rf=[R ; Ra];

%% Forcing terms
% Should be written as if the forcing angular frequency value is 1
% i.e. the forcing period is 2*pi
Forcing = zeros(2*sys.H+1,sys.nz_tot); % DO NOT CHANGE this line.

% if the equation number k is forced, write the forcing in Forced(:,k) as :
%Forcing(:,k) = forcing_function(t); 

end