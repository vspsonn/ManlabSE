function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% branches of solution.


if sys.H == 0
    x = Section.Upp(1,:);
    y = Section.Upp(2,:);
    lambda = Section.Upp(3,:);
    
    eig_val = Section.Eigen_init.values; % eigen values at the beginning of the solution-branch.
    
    figure(6)
    subplot(2,2,1);hold on;
    plot(lambda(1),real(eig_val),'bo');
    xlabel('\lambda');
    ylabel('Real part of the eigenvalues');
    
    subplot(2,2,3);hold on;
    plot(lambda(1),imag(eig_val),'bo');
    xlabel('\lambda');
    ylabel('Imaginary part of the eigenvalues');
    
    subplot(2,2,[2 4]);hold on;
    plot(real(eig_val),imag(eig_val),'bo');
    xlabel('Real part of the eigenvalues');
    ylabel('Imaginary part of the eigenvalues');
    
else
    
    figure(5);hold on;
    plotdiagnormHBM(sys,Section,1,'lambda');
    
    
end

end

