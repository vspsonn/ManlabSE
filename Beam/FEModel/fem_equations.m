function [Rf, dRf] = fem_equations(sys, Uf, dUf)


[R, Ra] = sys.parameters.assemble(Uf);
Rf =[R; Ra ];

dR    = zeros(sys.neq,1);
dRa   = zeros(sys.neq_aux,1);
dRf = [dR; dRa];

end