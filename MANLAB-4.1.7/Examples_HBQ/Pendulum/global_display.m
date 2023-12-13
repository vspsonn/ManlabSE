function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% User-defined function to display the solution-branch.
%
% Section is a 'CheckPoint class' object 
% containing all the information about the series:
% discrete representation, bifurcation analysis,...
%
% See the repertory MANLAB-4.x/SRC/Display for more functions.

% Upp = Section.Upp;              % Point representation of the series
% H=sys.H;                        % Number of harmonics
% lambda_sec = Upp(sys.getcoord('lambda'),:);  % Values of lambda along the Section
% omega_sec = Upp(sys.getcoord('omega'),:);    % Values of omega  along the Section


figure(5);hold on;

% Plot the L2 norm of the variables 1,2,3 with respect to omega.
plotbranchnormHBM(sys,Section,[1 2 3],'omega');

%%% Plot the norm of the first harmonic of the two first variables with
%%% respect to omega
%plotbranchHBM(sys,Section,[1 2],[1 1],'omega');

%%% Plot the norm and the phase of the first harmonic of the two first 
%%% variables with respect to omega
% plotbranchampphaseHBM(sys,Section,[1 2],[1 1],'omega');

%%% Plot the Floquet multipliers with respect to 'omega'
% plotfloquetmultHBM(sys,Section,'omega');

%%% Plot the Floquet exponents with respect to 'omega'
% plotfloquetexpHBM(sys,Section,'omega');

%%% Compute the first and the second sine and cosine coefficients of the
%%% first variable and put then in the matrices Ampcos/Ampsin.
% [Ampcos,Ampsin] = calcbranchsincosHBM(sys,Section,[1 1],[1 2]);

%%% Compute the phase of the first and the second harmonics of the first
%%% variable and put it in the matrix Phase.
% Phase = calcbranchphaseHBM(sys,Section,[1 1],[1 2]);

%%% Compute the amplitude of the first and the second harmonics of the
%%% first variable and put it in the matrix Hampl.
% Hampl = calcbranchHBM(sys,Section,[1 1],[1 2]);

%%% Compute the maximum and the minimum of the first and the second
%%% variable over the Section and put it in Max/Min matrices.
% [Max,Min]=calcmaxminHBM(sys,Section,[1 2])

end
