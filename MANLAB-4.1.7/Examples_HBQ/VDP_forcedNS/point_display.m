function [] = point_display(sys,U)

figure(11);clf;hold on;
% Plot the two first variables in the time domain over 1 period.
plotperiodHBM(sys,U,[1 2]);

figure(12);clf;hold on;
% Phase diagram (z1,z2) = (u,v).
plotphasediagHBM(sys,U,[1 2]);

end
