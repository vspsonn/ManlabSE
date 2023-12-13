function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% User-defined function to display solutions along the branch.
%
% See the repertory MANLAB-4.x/SRC/Display for more informations.

figure(11);clf; hold on;
% Plot the variable 1,2 over one period.
plotperiodHBM(sys,Uf,[1 2]);

figure(12);clf; hold on;
% Plot the variable 1 with respect to variable 2 
% over one period. [phase diagram] 
plotphasediagHBM(sys,Uf,[1 2]);

%%% Plot the amplitude of the sine and the cosine coefficients of the
%%% first and the second variable.
% plotbarsincosHBM(sys,Uf,[1 2]);

%%% Plot the amplitude of the harmonics of the first and the second
%%% variable.
% plotbarHBM(sys,Uf,[1 2])

%%% Compute the evolution in the time domain of the first and the second
%%% variable over 10 periods (the time is dimensionless) with 1e4 points.
% nb_period = 10;
% nb_samples = 1e4;
% time = linspace(0,nb_period,nb_samples);
% Utime = calcperiodHBM(sys,Uf,[1 2],time)

end
