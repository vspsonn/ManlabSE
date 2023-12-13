function [Res] = R(obj,U) 
% Compute the residu of the system at point U.

if isa(U,'Taylor')
 Res=Taylor(get(U,'order'),zeros(obj.neq,1));
else
 Res=zeros(obj.neq,1);
end

%% Equations defined by the user
u1=U(1); u2=U(2); lambda=U(3);

Res(1) = 2*u1 - u2 + 100*u1/(1+u1+u1*u1) - lambda;
Res(2) = 2*u2 - u1 + 100*u2/(1+u2+u2*u2) - (lambda+obj.parameters.mu);
  
