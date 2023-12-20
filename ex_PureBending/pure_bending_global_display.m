function [] = pure_bending_global_display(sys,Section)

fe_model = sys.parameters;

beam_def = fe_model.list_components{1};
ref_node = ComponentPartDef(beam_def, -1);

tracked_dofs = fe_model.get_nodal_dofs(ref_node, 1);
x = Section.Upp(tracked_dofs(3:4), :);

tracked_dofs = fe_model.get_nodal_dofs(ref_node, 0);
x0 = fe_model.U0(tracked_dofs(3:4));

lambda = Section.Upp(sys.neq+1,:);

x_analytical = [sin(x0(1) * lambda)./lambda; (1. - cos(x0(1) * lambda))./lambda];

figure(5)
plot(lambda, x(1, :) - x0(1), '-ob',...
    lambda, x(2, :) - x0(2), '-or');
hold on;
plot(lambda, x_analytical(1, :) - x0(1), '+b',...
    lambda, x_analytical(2, :) - x0(2), '+r');
legend('x', 'y', 'x_{anal}', 'y_{anal}', 'Location','NorthEastOutside');
xlabel('lambda');

end

