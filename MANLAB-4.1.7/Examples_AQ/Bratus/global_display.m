function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% branches of solution.

x1 = Section.Upp(1,:);
xend = Section.Upp(floor(sys.neq/2),:);
lambda = Section.Upp(sys.neq+1,:);

figure(5)
plot(lambda,x1,'-ob',lambda,xend,'-xr');hold on;
legend('x_1','x_{N/2}'); xlabel('lambda');

end

