classdef ComponentDef < handle
    properties
        n_dim
        number_of_elements
        list_elements
        list_node_numbers
    end
    
    methods
        function obj = ComponentDef(n_dim)
            obj.n_dim = n_dim;
            obj.list_elements = cell(1, 0);
            obj.list_node_numbers = [];
        end
        
        function n_dim = get_n_dim(obj)
            n_dim = obj.n_dim;
        end
        
        function [obj, U, ind_nodes] = mesh(obj, ind_nodes)
            U = [];
        end
    end
    
end