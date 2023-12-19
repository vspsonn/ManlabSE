classdef UnitQuatConstraint < Element
    properties
        n_dim
    end
    methods
        function obj = UnitQuatConstraint(node_number, n_dim)
            obj.n_dim = n_dim;
            obj.number_of_nodes = 1;
            
            if n_dim == 2
                obj.loc_dof = 1 + 4 * (node_number-1);
                obj.loc_config = (1:2) + 4 * (node_number-1);
            else
                obj.loc_dof = 1 + 7 * (node_number-1);
                obj.loc_config = (1:4) + 7 * (node_number-1);
            end
        end
        
        function [R, Ra] = assemble(obj, U, ~)
            n_cols = size(U, 2);
            
            Ra = zeros(0, n_cols);
            
            if obj.n_dim == 2
                R = (U(1, :) .* U(1, :) + U(2, :) .* U(2, :)) - ones(1, n_cols);
            else
                R = (U(1, :) .* U(1, :) + U(2, :) .* U(2, :) + U(3, :) .* U(3, :) + U(4, :) .* U(4, :)) - ones(1, n_cols);
            end
        end
    end
end

