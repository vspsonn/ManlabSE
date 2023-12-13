classdef UnitQuatConstraint < Element

    properties
        n_dof
    end

    methods
        function obj = UnitQuatConstraint(node_number, n_dof)
            obj.n_nodes = 1;
            obj.loc_dof = 1 + (n_dof+1) * (node_number-1);
            obj.loc_config = (1:4) + (n_dof+1) * (node_number-1);

            obj.n_dof = n_dof;
        end

        function [R, Ra] = assemble(obj, U, ~)
            Ra = zeros(0, 1);

            R = U(1:4)' * U(1:4) - 1.0;
        end
    end
end

