classdef ComponentDef
    properties
        n_elem
        list_elem
    end

    methods
        function [obj, U, ind_nodes] = mesh(obj, ind_nodes)
            U = [];
        end
    end

end