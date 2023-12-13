function [Rf,dRa,Forced] = equations(sys,t,zf,dzf,d2zf)

% Pendulum : x'' + lambda x' + sin(x) = 0
% x'' + lambda x' + s = 0
% s' =       x' * c
% c' =     - x' * s
% U = [x,lambda,s,c];
%%% s(0) = sin(x(0))
%%% c(0) = cos(x(0))

x = zf(1); dx = dzf(1); d2x = d2zf(1);
s = zf(2); ds = dzf(2);
c = zf(3); dc = dzf(3);
lambda = zf(end);

%%% Main part
R = d2x + lambda*dx + s;   

%%% Auxiliary part
Ra = zeros(sys.nz_aux,1);     dRa = zeros(sys.nz_aux,1);
Ra(1) = s - sin(x);           dRa(1) = ds - dx*c;
Ra(2) = c - cos(x);           dRa(2) = dc + dx*s;     

%%% Concatenation
Rf = [R;Ra];

%%% Forced terms :
Forced = zeros(2*sys.H+1,sys.nz_tot);

end

