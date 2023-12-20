classdef ExternalForceDef < ComponentDef
    properties
        component_part
        force_props
    end

    methods
        function obj = ExternalForceDef(component_part, force_props)
            obj = obj@ComponentDef(component_part.get_n_dim());
            obj.number_of_elements = 1;
            obj.component_part = component_part;
            
            obj.force_props = force_props;
        end
        
        function [obj, U, ind_nodes] = mesh(obj, ind_nodes)
            U = [];
            
            node_number = obj.component_part.get_node_number();
            obj.list_node_numbers = [node_number];
            obj.list_elements = {ExternalForceElement(node_number, obj.force_props, obj.n_dim)};
        end
    end
end