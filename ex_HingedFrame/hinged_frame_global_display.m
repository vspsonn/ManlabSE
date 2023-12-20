function [] = hinged_frame_global_display(sys,Section)

fe_model = sys.parameters;

beam2_def = fe_model.list_components{2};
ref_node = ComponentPartDef(beam2_def, 0);

tracked_dofs = fe_model.get_nodal_dofs(ref_node, 1);
x = Section.Upp(tracked_dofs(3:4), :);

tracked_dofs = fe_model.get_nodal_dofs(ref_node, 0);
x0 = fe_model.U0(tracked_dofs(3:4));

lambda = Section.Upp(sys.neq+1,:);

figure(5)
plot(lambda, x(1, :) - x0(1), '-r',...
    lambda, -(x(2, :) - x0(2)), '-b');
hold on;
legend('x', 'y'); 
xlabel('lambda');

end

