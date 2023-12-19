function [] = fem_point_display(sys, Uf)

figure
hold on
sys.parameters.plot(1, Uf);
axis equal
grid on

end

