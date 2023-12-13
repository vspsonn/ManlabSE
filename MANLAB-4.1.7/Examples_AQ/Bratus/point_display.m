function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% This is the function used to display during the computation of the
% branches of solution.

N=sys.parameters.N;        % number of mesh interval

x=[ 0 : 1/(N+1)  : 1 ];
u=[ 0 , Uf(1:N)' , 0] ;

figure(11)
plot(x,u)
title(['BRATU : u(x) for lambda = ' num2str(Uf(sys.neq+1))],'Fontsize',16)
xlabel('x');ylabel('u(x)');

end

