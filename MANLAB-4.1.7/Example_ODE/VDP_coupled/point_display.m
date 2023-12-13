function [] = point_display(sys,Uf)

H=sys.H;

switch sys.type
    case 'HBM'
        
        % Computation of the periodic solution over 1 period : time is a
        % dimensionless time vector.
        time = linspace(0,1,1000);
        Ut = calcperiodHBM(sys,Uf,[1 2 3 4],time);
        
    case 'QPHBM'
        % The quasi-periodic solution is reconstructed over the grid
        % time below from the Fourier coefficients.
        
        omega1 = Uf(sys.neq-1);
        omega2 = Uf(sys.neq);
        
        Tmax = pi/omega1 + pi/omega2;
        
        nb = 50;
        time=linspace(0,nb*Tmax,nb*400);
        Ut = calctimeQPHBM(sys,Uf,[1 2 3 4],time);
        
        % The matrix of coefficients of the double Fourier are displayed
        % for the first variable, q1, (top) and the third, q2, (bottom)
        figure(13);clf;
        subplot(2,1,1);
        plotcoefQPHBM(sys,Uf,1);
        subplot(2,1,2);
        plotcoefQPHBM(sys,Uf,3);
end

figure(11);
plot(time,Ut(:,1),'b',time,Ut(:,3),'r')
title(['Time evolution of q_1 and q_2 with ' num2str(H) ' harmonics'],'Fontsize',16)
xlabel('time (s)');legend(['q_1'],['q_2'],'Location','SouthEast');grid on;

figure(12)
plot3(Ut(:,1),Ut(:,2),Ut(:,3))
xlabel('q_1')
ylabel('p_1')
zlabel('q_2')

end
