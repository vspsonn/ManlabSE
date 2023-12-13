function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% This is the function used to display during the computation of the
% branches of solution.

u      = Uf(1:sys.neq);
lambda = Uf(sys.neq+1);
Ua     = Uf(sys.neq+2:end);

%%% Plot the value of the main variables and the value of the auxiliary
%%% variables
figure(11)
subplot(2,1,1)
bar(u)
title('Main variables');
subplot(2,1,2)
bar(Ua)
title('Auxiliary variables');

end

