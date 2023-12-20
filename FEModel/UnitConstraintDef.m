classdef UnitConstraintDef < ComponentDef
    methods
        function obj = UnitConstraintDef(n_dim, node_numbers)
            obj = obj@ComponentDef(n_dim);
            
            obj.list_node_numbers = node_numbers;
            obj.list_elements = {};
            for i = 1:size(node_numbers, 2)
                obj.list_elements{i} = UnitQuatConstraint(node_numbers(i), obj.n_dim);
            end
        end
    end
end

