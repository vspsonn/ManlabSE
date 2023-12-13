classdef BoundaryConditionDef < ComponentDef
    properties
        n_dof
        node_num
        node_LM
    end
    
    methods
        function obj = BoundaryConditionDef(node_num, n_dof)
            obj = obj@ComponentDef();
            obj.n_dof = n_dof;
            obj.n_elem = 1;
            
            obj.node_num = node_num;
        end
        
        function [obj, U, ind_nodes] = mesh(obj, ind_nodes)
            U = zeros(obj.n_dof,1);
            
            obj.node_LM = ind_nodes + 1;
            
            obj.list_elem = {BoundaryConditionElement(obj.node_num, obj.node_LM, obj.n_dof)};
            
            ind_nodes = ind_nodes + 1;
        end
    end
    
end