function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% solution branch.

% Section.U0  (Taylor series) contains the Taylor series used to describe the section.
% Section.Amax (real positive number) is the domain of validity of the series.
% Section.Upp (ninc,nb_pts) is the point representation of the series.

u1 = Section.Upp(1,:);
lambda = Section.Upp(sys.neq+1,:);

%%% plot the value of the first variable u1 with respect to the value of
%%% the continuation parameter lambda
figure(5)
plot(lambda,u1);hold on;

end

