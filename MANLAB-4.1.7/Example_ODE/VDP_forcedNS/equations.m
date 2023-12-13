function [Rf,dRa,Forcing] = equations(sys,t,zf,dzf,d2zf)
 
%    DAE for a forced van der pol model ( u' = u // v' - (mu1 - mu2 u - u^2) v + u = F cos(lambda t) )  (nz=2)
%
% The derivated terms should be with a -1 and the forced term is
% isolated on the right hand side, as follows :
%
%    v                    - u'  = 0
%  - u              + rv  - v'  = - F cos(lambda t)
%  mu1 - r - mu2*u  - u^2       = 0
%

%%% parameters of the system
mu1 = sys.parameters.mu1;
mu2 = sys.parameters.mu2;
F = sys.parameters.F;

%%% variables of the system
u = zf(1);du=dzf(1);
v = zf(2);dv=dzf(2);
r = zf(3);
lambda = zf(end);

%%% Physical equations
R =zeros(sys.nz,1);
R(1) = v - du;
R(2) = v*r - u - dv;

%%% Auxiliary equations
Ra = mu1 - r - mu2*u - u^2;

%%% All the equations of the quadratic system
Rf = [ R;Ra];
dRa = zeros(size(Ra));

%% Forced terms
% Should be written as if the forcing pulsation value is 1
Forcing = zeros(2*sys.H+1,sys.nz_tot);
Forcing(:,2) = -F*cos(t); % puls = 1 to write the forced terms.

end