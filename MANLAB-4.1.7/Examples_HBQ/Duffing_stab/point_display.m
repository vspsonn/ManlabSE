function [] = point_display(sys,U)

figure(11);clf;hold on;
% Plots the two variables in the time domain over 1 period.
plotperiodHBM(sys,U,[1 2]);

end
