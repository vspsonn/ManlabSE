function [] = global_display(sys,Section)
% function [] = global_display(sys,Utot)
% This is the function used to display during the computation of the
% branches of solution.

Upp = Section.Upp;

mu = Upp(sys.neq+1,:);
Xn = Upp(2:sys.neq,:);

figure(5);hold on;
plot(mu,Xn(end-1,:),'b');
plot(mu,Xn(end,:),'r');

end

