classdef FEModel
    
    properties
        n_dim
        
        list_components
        list_components_bc_rot
        list_components_bc_trans
        
        list_config
        list_dof_free
        list_dof_fix
        
        n_dof
        n_dof_free
        n_dof_fix
        U0
        
        n_aux
        Ua0
    end
    
    methods
        function obj = FEModel(n_dim)
            obj.n_dim = n_dim;
            
            obj.list_components = cell(1,0);
            obj.list_components_bc_rot = cell(1,0);
            obj.list_components_bc_trans = cell(1,0);
        end
        
        function obj = add_component(obj, component)
            obj.list_components{size(obj.list_components, 2)+1} = component;
        end
        
        function obj = add_bc_rot(obj, component_part)
            obj.list_components_bc_rot{size(obj.list_components_bc_rot, 2)+1} = component_part;
        end
        
        function obj = add_bc_trans(obj, component_part)
            obj.list_components_bc_trans{size(obj.list_components_bc_trans, 2)+1} = component_part;
        end
        
        function obj = process_bc(obj, ind_nodes)
            obj.n_dof = size(obj.U0, 1);
            
            obj.list_dof_free = [];
            obj.list_dof_fix = [];
            nodes_fix_rot = [];
            if obj.n_dim == 2
                for i = 1:size(obj.list_components_bc_trans, 2)
                    num_node = obj.list_components_bc_trans{i}.get_node_number();
                    obj.list_dof_fix = [obj.list_dof_fix (3:4) + 4*(num_node-1)];
                end
                for i = 1:size(obj.list_components_bc_rot, 2)
                    num_node = obj.list_components_bc_rot{i}.get_node_number();
                    obj.list_dof_fix = [obj.list_dof_fix (1:2) + 4*(num_node-1)];
                    nodes_fix_rot = [nodes_fix_rot num_node];
                end
            else
                for i = 1:size(obj.list_components_bc_trans, 2)
                    num_node = obj.list_components_bc_trans{i}.get_node_number();
                    obj.list_dof_fix = [obj.list_dof_fix (5:7) + 7*(num_node-1)];
                end
                for i = 1:size(obj.list_components_bc_rot, 2)
                    num_node = obj.list_components_bc_rot{i}.get_node_number();
                    obj.list_dof_fix = [obj.list_dof_fix (1:4) + 7*(num_node-1)];
                    nodes_fix_rot = [nodes_fix_rot num_node];
                end
            end
            
            nodes_free_rot = zeros(1, ind_nodes - size(nodes_fix_rot, 2));
            ind_free = 1;
            for i = 1:ind_nodes
                ind = find(nodes_fix_rot == i, 1);
                if isempty(ind)
                    nodes_free_rot(ind_free) = i;
                    ind_free = ind_free + 1;
                end
            end
            obj.list_components{size(obj.list_components, 2)+1} = UnitConstraintDef(obj.n_dim, nodes_free_rot);
            
            obj.n_dof_fix = size(obj.list_dof_fix, 2);
            obj.n_dof_free = obj.n_dof - obj.n_dof_fix;
            obj.list_dof_free = zeros(1, obj.n_dof_free);
            ind_free = 1;
            for i = 1:obj.n_dof
                ind = find(obj.list_dof_fix == i, 1);
                if isempty(ind)
                    obj.list_dof_free(ind_free) = i;
                    ind_free = ind_free + 1;
                end
            end
        end
        
        function obj = mesh(obj)
            obj.U0 = [];
            ind_nodes = 0;
            for i = 1:size(obj.list_components, 2)
                [obj.list_components{i}, U0i, ind_nodes] = obj.list_components{i}.mesh(ind_nodes);
                obj.U0 = [obj.U0; U0i];
            end
            
            obj = obj.process_bc(ind_nodes);
        end
        
        function obj = initialize(obj)
            obj.Ua0 = [];
            ind_aux = 1;
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elements, 2)
                    elem = obj.list_components{i}.list_elements{j};
                    [obj.list_components{i}.list_elements{j}, Ua0_j, ind_aux] = elem.initialize(obj, obj.U0(elem.loc_config), ind_aux);
                    obj.Ua0 = [obj.Ua0; Ua0_j];
                end
            end
            
            obj.n_aux = size(obj.Ua0, 1);
        end
        
        function u0 = get_U0(obj, continuation_parameter)
            u0 = [obj.U0(obj.list_dof_free); continuation_parameter; obj.Ua0];
        end
        
        function dofs = get_nodal_dofs(obj, component_part, in_free)
            node_number = component_part.get_node_number();
            
            if obj.n_dim == 2
                dofs = [(1:4) + 4*(node_number-1)];
            else
                dofs = [(1:7) + 7*(node_number-1)];
            end
            
            if in_free
                full_dofs = dofs;
                for i=1:size(full_dofs, 2)
                    ind = find(obj.list_dof_free == full_dofs(i), 1);
                    if isempty(ind)
                        dofs(i) = 0;
                    else
                        dofs(i) = ind;
                    end
                end
            end
        end
        
        function [R, Ra] = assemble(obj, Uf)
            n_col = size(Uf, 2);
            U = zeros(obj.n_dof+1, n_col);
            for i = 1:n_col
                U(1:obj.n_dof, i) = obj.U0;
            end
            U(obj.list_dof_free, :) = Uf(1:obj.n_dof_free, :);
            U(obj.n_dof+1, :) = Uf(obj.n_dof_free+1, :);
            Ua = Uf(obj.n_dof_free+2:end, :);
            
            R = zeros(obj.n_dof_free, n_col);
            Ra = zeros(obj.n_aux, n_col);
            
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elements, 2)
                    elem = obj.list_components{i}.list_elements{j};
                    loc_aux = elem.loc_aux + (0:(elem.number_of_aux-1));
                    
                    [R_i, Ra_i] = elem.assemble(U(elem.loc_config, :), Ua(loc_aux, :));
                    
                    loc_dof_free = elem.get_loc_dof_free(obj.list_dof_free);
                    for m = 1:size(loc_dof_free, 2)
                        if loc_dof_free(m) ~= 0
                            R(loc_dof_free(m), :) = R(loc_dof_free(m), :) + R_i(m, :);
                        end
                    end
                    Ra(loc_aux, :) = Ra_i;
                end
            end
        end
        
        function plot(obj, fig, Uf)
            if fig == 0
                figure
                hold on
            end
            
            U = [obj.U0; 0];
            U(obj.list_dof_free) = Uf(1:obj.n_dof_free);
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elements, 2)
                    elem = obj.list_components{i}.list_elements{j};
                    elem.plot(fig, U(elem.loc_config));
                end
            end
            
            if fig == 0
                axis equal
                grid on
            end
        end
        
    end
end