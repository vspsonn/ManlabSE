function [Rf,dRa,Forcing] = equations(sys,t,zf,dzf,d2zf)
 
% The derivated terms should be with a -1 sign and the forced term is
% isolated on the right hand side, as follows :
%
% p1' = q1
% q1' = - w1^2*p1-a1*q1+a2*lambda*(q1+q2)-a3*lambda*w-a4*lambda*v 
% p2' = q2
% q2' =- w1^2*p2-b1*q2+b2*lambda*(q1+q2)-b3*lambda*w-b4*lambda*v
%  0   =r - (p1 + p2)^2 
%  0   =w - (p1+p2)*(q1+q2)            
%  0   =v- r*(q1+q2) 

%%% parameters of the system
w1 = sys.parameters.w1;
w2 = sys.parameters.w2;

a1 = sys.parameters.a1;
a2 = sys.parameters.a2;
a3 = sys.parameters.a3;
a4 = sys.parameters.a4;

b1 = sys.parameters.b1;
b2 = sys.parameters.b2;
b3 = sys.parameters.b3;
b4 = sys.parameters.b4;

%%% variables of the system
p1 = zf(1); dp1=dzf(1);
q1 = zf(2); dq1=dzf(2);
p2 = zf(3); dp2=dzf(3);
q2 = zf(4); dq2=dzf(4);
r = zf(5);
w = zf(6);
v = zf(7);
lambda = zf(end);

%%% Physical equations
R = zeros(sys.nz,1);      
R(1) = q1 - dp1;
R(2) = -w1^2*p1-a1*q1+lambda*(a2*(q1+q2)-a3*w-a4*v) - dq1;
R(3) = q2 - dp2;
R(4) = -w2^2*p2-b1*q2+lambda*(b2*(q1+q2)-b3*w-b4*v) - dq2;

%%% Auxiliary equations
Ra = zeros(sys.nz_aux,1); dRa = zeros(sys.nz_aux,1);
Ra(1) = r - (p1 + p2)^2;
Ra(2) = w - (p1+p2)*(q1+q2);
Ra(3) = v - r*(q1+q2);

%%% All the equations of the quadratic system
Rf = [ R;Ra];        

%%% Forcing terms
Forcing = zeros(2*sys.H+1,sys.nz_tot);

end