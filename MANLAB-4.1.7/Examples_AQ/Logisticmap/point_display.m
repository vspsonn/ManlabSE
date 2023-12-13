function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% This is the function used to display during the computation of the
% branches of solution.

x = Uf(1:sys.neq);
mu= Uf(sys.neq+1);

nb = min(30,sys.neq);

figure(11);

subplot(3,1,1);
plot(x(1:nb),'-o');
title([num2str(nb) ' first points of the sequence for mu = ' num2str(mu)]);
ax = gca;
set(ax,'Fontsize',18);

subplot(3,1,2);
plot(x(end-nb+1:end),'-x');
ylim([0 1]);
title([num2str(nb) ' last points of the sequence for mu = ' num2str(mu)]);
ax = gca;
set(ax,'Fontsize',18);

subplot(3,1,3);
plot(x,'.');
title([num2str(sys.neq) ' first points of the sequence for mu = ' num2str(mu)]);
ax = gca;
set(ax,'Fontsize',18);

end

