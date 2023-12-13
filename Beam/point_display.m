function [] = point_display(sys, Uf)

figure
hold on
sys.parameters.plot([], Uf);
axis equal
grid on

end

