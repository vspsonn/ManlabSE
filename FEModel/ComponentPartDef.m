classdef ComponentPartDef
    properties
        base
        node_index_on_base
    end
    
    methods
        function obj = ComponentPartDef(base, node_index_on_base)
            obj.base = base;
            obj.node_index_on_base = node_index_on_base;
        end
        
        function n_dim = get_n_dim(obj)
            n_dim = obj.base.get_n_dim();
        end
        
        function node_number = get_node_number(obj)
            if obj.node_index_on_base == -1
                node_number = obj.base.list_node_numbers(end);
            else
                if obj.node_index_on_base < 1
                    index = floor(obj.node_index_on_base * (size(obj.base.list_node_numbers, 2)-1) + 1);
                    node_number = obj.base.list_node_numbers(index);
                else
                    node_number = obj.base.list_node_numbers(obj.node_index_on_base);
                end
            end
        end
    end
end

