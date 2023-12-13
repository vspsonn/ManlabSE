function [] = global_display(sys,Section)
% Section is an object
% containing all the information about the series:
% discrete representation, bifurcation analysis,...


% Amplitude and phase of the first and the second harmonics 
% of the first variable with respect to omega in figure 5
figure(5);hold on;
plotbranchampphaseHBM(sys,Section,[1 1],[1 3],'omega');

% Plots Floquet exponents/multipliers.
% figure(7);hold on;
% plotfloquetexpHBM(sys,Section,'omega');

figure(8);hold on;
plotfloquetmultHBM(sys,Section,'omega');


end
