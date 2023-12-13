function [] = point_display(sys,Utot)
% function [] = point_display(sys,Utot)
% This is the function used to display during the computation of the
% branches of solution.

u      = Utot(1:sys.neq);
%lambda = Utot(sys.neq+1);
%Uaux   = Utot(sys.neq+2:end);

figure(11);
bar(u)

end

