function [] = bend45_global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% branches of solution.

fe_model = sys.parameters;

beam_def = fe_model.list_components{1};
ef_node = ComponentPartDef(beam_def, -1);

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 1);
x = Section.Upp(tracked_dofs(5:7), :);

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 0);
x0 = fe_model.U0(tracked_dofs(5:7));

lambda = Section.Upp(sys.neq+1,:);

figure(5)
plot(lambda, x(1, :) - x0(1), '-ob',...
    lambda, x(2, :) - x0(2), '-xr',...
    lambda, x(3, :) - x0(3), '-+g');
hold on;
legend('x','y','z');
xlabel('lambda');

end

