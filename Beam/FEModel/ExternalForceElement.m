classdef ExternalForceElement < Element
    
    properties
        n_dim
        n_dof
        force_props
    end
    
    methods
        function obj = ExternalForceElement(node_number, force_props, n_dim)
            obj = obj@Element();
            
            obj.n_dim = n_dim;
            obj.number_of_nodes = 1;
            
            if n_dim == 2
                obj.n_dof = 3;
                obj.loc_dof = [3:4 2] + 4 * (node_number-1);
                
                n_config = 4;
            else
                obj.n_dof = 6;
                obj.loc_dof = [5:7 2:4] + 7 * (node_number-1);
                
                n_config = 7;
            end
            
            obj.loc_config = zeros(1, n_config);
            obj.loc_config = (1:n_config) + n_config * (node_number-1);
            
            obj.force_props = force_props;
        end
        
        function [obj, Ua, loc_aux] = initialize(obj, fe_model, ~, loc_aux)
            obj.loc_config = [obj.loc_config fe_model.n_dof+1];
            
            if obj.force_props.follower
                obj.number_of_aux = 0;
                obj.loc_aux = [];
                Ua = [];
            else
                if obj.n_dim == 2
                    obj.number_of_aux = 2;
                else
                    obj.number_of_aux = 6;
                end
                
                obj.loc_aux = loc_aux;
                loc_aux = loc_aux + obj.number_of_aux;
                Ua = zeros(obj.number_of_aux, 1);
            end
        end
        
        function [R, Ra] = assemble(obj, U, Ua)
            n_cols = size(U, 2);
            
            n_config = obj.n_dof + 1;
            F = - obj.force_props.dir * U(n_config+1, :);
            
            if obj.force_props.follower
                Ra = zeros(0, n_cols);
                R = F;
                
            else
                R = zeros(obj.n_dof, n_cols);
                
                if obj.n_dim == 2
                    q = U(2, :);
                    
                    yu = Ua(1:2, :);
                    Ra = yu - q .* [-F(2, :); F(1, :)];
                    
                    q0 = U(1, :);
                    R(1:2, :) = F(1:2, :) - (2.0 * q0) .* yu + (2.0 * q) .* [-yu(2, :); yu(1, :)];
                    R(3, :) = F(3, :);
                else
                    q = U(2:4, :);
                    
                    Ra = zeros(6, n_cols);
                    yu = Ua(1:3, :);
                    Ra(1:3, :) = yu - SO3_cross(q, F(1:3, :));
                    yr = Ua(4:6, :);
                    Ra(4:6, :) = yr - SO3_cross(q, F(4:6, :));
                    
                    q0 = U(1, :);
                    R(1:3, :) = F(1:3, :) - 2.0 * q0 .* yu + 2.0 * SO3_cross(q, yu);
                    R(4:6, :) = F(4:6, :) - 2.0 * q0 .* yr + 2.0 * SO3_cross(q, yr);
                end
            end
        end
    end
end
