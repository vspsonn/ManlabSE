classdef BoundaryConditionElement < Element
    properties
        n_dof
    end

    methods
        function obj = BoundaryConditionElement(node_num, node_LM, n_dof)
            obj.n_nodes = 1;
            obj.n_dof = n_dof;
            n_config = n_dof+1;
            obj.loc_config = [(1:n_config) + n_config * (node_num-1) (1:n_dof)+(n_dof+1)*(node_LM-1)];
            obj.loc_dof = [(2:(n_dof+1)) + (n_dof+1) * (node_num-1) (1:n_dof)+(n_dof+1)*(node_LM-1)];
        end

        function [R, Ra] = assemble(obj, U, Ua)
            Ra = zeros(0, 1);

            R = zeros(2 * obj.n_dof, 1);
            R(1:obj.n_dof) = U(obj.n_dof+1 + (1:obj.n_dof));
            R(obj.n_dof + (1:obj.n_dof)) = U(1 + (1:obj.n_dof));
        end
    end
end

