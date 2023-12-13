function [] = point_display(sys,Uf)

H=sys.H;

switch sys.type
    case 'HBM'
        
        % Computation of the periodic solution over 1 period : time is a
        % dimensionless time vector.
        time = linspace(0,1,1000);
        Ut = calcperiodHBM(sys,Uf,[1 2],time);
        
    case 'QPHBM'
        % The quasi-periodic solution is reconstructed over the grid
        % time below from the Fourier coefficients.
        
        omega1 = Uf(sys.neq-1);
        omega2 = Uf(sys.neq);
        
        Tmax = pi/omega1 + pi/omega2;
        
        nb = 50;
        time=linspace(0,nb*Tmax,nb*400);
        Ut = calctimeQPHBM(sys,Uf,[1 2],time);
        
        % The matrix of coefficients of the double Fourier are displayed
        % for the first variable, u, (top) and the second, v, (bottom)
        figure(13);clf;
        subplot(2,1,1);
        plotcoefQPHBM(sys,Uf,1);
        subplot(2,1,2);
        plotcoefQPHBM(sys,Uf,2);
end

figure(11);
plot(time,Ut(:,1),'b',time,Ut(:,2),'r')
title(['Time evolution of u and v with ' num2str(H) ' harmonics'],'Fontsize',16)
xlabel('time (s)');legend(['u'],['v'],'Location','SouthEast');grid on;

figure(12)
plot(Ut(:,1),Ut(:,2))
xlabel('u')
ylabel('v')

end
