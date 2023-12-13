classdef Element
    properties
        n_nodes
        loc_config
        loc_dof
        
        n_aux
        loc_aux
    end
    
    methods
        function [obj, Ua, loc_aux] = initialize(obj, fe_model, U, loc_aux)
            obj.n_aux = 0;
            obj.loc_aux = [];
            Ua = [];
        end
        
        function plot(obj, fig, U)
        end
    end
end

