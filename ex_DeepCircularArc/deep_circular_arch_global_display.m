function [] = deep_circular_arch_global_display(sys,Section)
% function [] = global_display(sys,Section)
% This is the function used to display during the computation of the
% branches of solution.

fe_model = sys.parameters;

beam_def = fe_model.list_components{1};
ef_node = ComponentPartDef(beam_def, 0.5);

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 1);
x = Section.Upp(tracked_dofs(4), :);

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 0);
x0 = fe_model.U0(tracked_dofs(4));

lambda = Section.Upp(sys.neq+1,:);

figure(5)
plot(lambda, -(x - x0), '-b');
hold on;
legend('x'); 
xlabel('lambda');

end

