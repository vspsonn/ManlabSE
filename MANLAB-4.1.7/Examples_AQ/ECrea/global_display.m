function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% branches of solution.

u1 = Section.Upp(1,:);
u2 = Section.Upp(2,:);
lambda = Section.Upp(3,:);

figure(5)
plot(lambda,u1,'b',lambda,u2,'r');hold on;

end

