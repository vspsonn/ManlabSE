function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% This is the function used to display during the computation of the
% branches of solution.

u      = Uf(1:sys.neq);
lambda = Uf(sys.neq+1);
Ua   = Uf(sys.neq+2:end);

figure(11)
bar(u)

end

