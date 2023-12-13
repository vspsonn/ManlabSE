function [] = global_display(sys,Section)
% Section is a 'CheckPoint class' object 
% containing all the information about the series:
% discrete representation, bifurcation analysis,...

% Upp = Section.Upp; % Point representation of the series
% H=sys.H;
% nz_tot=sys.nz_tot;

figure(6);hold on;
% L2 norm of the first variable with respect to lambda in figure 6
plotbranchnormHBMbif(sys,Section,1,'lambda');

% Plot Floquet exponents and multipliers.
figure(7);hold on;
plotfloquetexpHBM(sys,Section,'lambda');

%figure(8);hold on;
%plotfloquetmultHBM(sys,Section,'lambda');


end
