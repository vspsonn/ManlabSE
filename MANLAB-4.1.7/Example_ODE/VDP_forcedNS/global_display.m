function [] = global_display(sys,Sec_diag)
% Section is a structure containing all the information of the series,
% its change of stability, simple bifurcation information,...

% For now, no functions has been written for the display of
% quasi-periodic branches of solutions.

Diag = calcdiagUpp(sys,Sec_diag); % Compute a structure containing the points over the whole diagram

switch sys.type
    case 'HBM'
        % Plot the L2 norm of the u (first variable) and v (second variable)
        % along the solution-branch with respect to lambda.
        figure(5);hold on;
        plotdiagnormHBMbif(sys,Diag,[1 2],'lambda');
        
        % Plot the Floquet exponents along the solution-branch.
        figure(7);hold on;
        plotfloquetexpHBM(sys,Sec_diag,'lambda');
        
    case 'QPHBM'
        % Plot the L2 norm of the u (first variable) and v (second variable)
        % along the solution-branch with respect to lambda.
        figure(5);hold on;
        plotdiagnormQPHBM(sys,Diag,[1 2],'lambda');
        
        omega1 = Diag.DiagUpp(sys.neq-1,:);
        omega2 = Diag.DiagUpp(sys.neq,:);
        lambda = Diag.DiagUpp(sys.neq+1,:);

        figure(6);hold on;
        plotdiagXYQPHBM(sys,Diag,lambda,omega2./omega1);
        xlabel('\lambda')
        ylabel('\omega_2 / \omega_1')
end

end
