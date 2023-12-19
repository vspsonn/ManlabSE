classdef Element
    properties
        number_of_nodes
        loc_config
        loc_dof
        
        number_of_aux
        loc_aux
    end
    
    methods
        function [obj, Ua, loc_aux] = initialize(obj, fe_model, U, loc_aux)
            obj.number_of_aux = 0;
            obj.loc_aux = [];
            Ua = [];
        end
        
        function loc_dof_free = get_loc_dof_free(obj, FEModel_loc_dof_free)
            loc_dof_free = zeros(1, size(obj.loc_dof, 2));
            for i = 1:size(obj.loc_dof, 2)
                ind = find(FEModel_loc_dof_free == obj.loc_dof(i), 1);
                if ~isempty(ind)
                    loc_dof_free(i) = ind;
                end                
            end            
        end
        
        function plot(obj, fig, U)
        end
    end
end

