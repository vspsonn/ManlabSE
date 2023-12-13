function [Rf,dRa,Forcing] = equations(sys,t,zf,dzf,d2zf)
 
%    DAE for the forced duffing model ( u' = v // v' = - epsilon*v - u - u^3 = alpha*cos(lambda t)  (nz=2)
%
% The derivated terms should be with a -1 sign and the forced term is
% isolated on the right hand side, as follows :
%
%    v                   - u'  = 0
%  - u - epsilon v - ru  - v'  = - alpha cos(lambda t)
%    r             - u^2       = 0

%%% parameters of the system
alpha = sys.parameters.alpha;
epsilon = sys.parameters.epsilon;

%%% variables of the system
u = zf(1);du=dzf(1);
v = zf(2);dv=dzf(2);
r = zf(3);
lambda = zf(end);

%%% Physical equations
R = zeros(sys.nz,1);         
R(1) = v - du;
R(2) = - u - epsilon*v - r*u - dv;

%%% Auxiliary equations
Ra = zeros(sys.nz_aux,1);     dRa = zeros(sys.nz_aux,1);
Ra(1) = r - u^2;

%%% All the equations of the quadratic system
Rf = [R;Ra];               

%%% Forced terms
% Should be written as if the forcing pulsation value is 1
% i.e. the forcing period is 2*pi
Forcing = zeros(2*sys.H+1,sys.nz_tot);
Forcing(:,2) = -alpha*cos(t);

end