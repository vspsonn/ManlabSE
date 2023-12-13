classdef ExternalForceDef < ComponentDef

    properties
        node_number
        force_props
    end

    methods
        function obj = ExternalForceDef(node_number, force_props)
            obj = obj@ComponentDef();

            obj.node_number = node_number;
            obj.force_props = force_props;
        end

        function [obj, U_def, ind_nodes] = mesh(obj, ind_nodes)
            obj.list_elem = {ExternalForceElement(obj.node_number, obj.force_props)};
            U_def = [];
        end
    end

end