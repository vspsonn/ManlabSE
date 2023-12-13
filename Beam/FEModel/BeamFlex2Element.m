classdef BeamFlex2Element < Element
    
    properties
        beam_props
        
        n_gp = 2
        J0
        F0
    end
    
    properties (Constant)
        loc_gp = {[0.], [-1/sqrt(3), 1/sqrt(3)], [0.0, -sqrt(3/5), sqrt(3/5)]};
        weight_gp = {[2.], [1., 1.], [8.0/9.0, 5.0/9.0, 5.0/9.0]};
        
        shape_function = {[],
            @(eta) [0.5*(1-eta), 0.5*(1+eta)],
            @(eta) [0.5*eta*(eta-1), 1-eta*eta, 0.5*eta*(eta+1)],
            @(eta) [-9/16*(1-eta)*(1/3+eta)*(1/3-eta), 27/16*(1-eta)*(1+eta)*(1/3-eta), 27/16*(1-eta)*(1+eta)*(1/3+eta), -9/16*(1+eta)*(1/3+eta)*(1/3-eta)]};
        der_shape_function = {[],
            @(eta) [-0.5, 0.5],
            @(eta) [eta-0.5, -2*eta, eta+0.5],
            @(eta) [-1.6875*eta*eta+1.125*eta+0.0625, 5.0625*eta*eta-1.125*eta-1.6875, -5.0625*eta*eta-1.125*eta+1.6875, 1.6875*eta*eta+1.125*eta-0.0625]};
    end
    
    methods
        function obj = BeamFlex2Element(list_node_numbers, beam_props)
            
            obj.n_nodes = size(list_node_numbers, 2);
            
            obj.loc_dof = zeros(1, 3 * obj.n_nodes);
            obj.loc_config = zeros(1, 4 * obj.n_nodes);
            for i = 1:obj.n_nodes
                obj.loc_dof((1:3) + (i-1) * 3) = (2:4) + 4 * (list_node_numbers(i)-1);
                obj.loc_config((1:4) + (i-1) * 4) = (1:4) + 4 * (list_node_numbers(i)-1);
            end
            
            obj.beam_props = beam_props;
        end
        
        function [obj, Ua, loc_aux] = initialize(obj, fe_model, U, loc_aux)
            
            obj.n_aux = obj.n_gp * (4 * obj.n_nodes + 3);
            obj.loc_aux = loc_aux;
            
            Ua = zeros(obj.n_aux, 1);
            obj.J0 = zeros(1, obj.n_gp);
            obj.F0 = zeros(3, obj.n_gp);
            
            ind_aux = 1;
            for i = 1:obj.n_gp
                eta = obj.loc_gp{obj.n_gp}(i);
                
                shpfcn = obj.shape_function{obj.n_nodes}(eta);
                der_shpfcn = obj.der_shape_function{obj.n_nodes}(eta);
                
                q_bar = zeros(4,1);
                for k = 1:obj.n_nodes
                    f_k = U((1:4) + (k-1)*4);
                    q_bar = q_bar + shpfcn(k) * f_k(1:4);
                end
                %z2 = q_bar' * q_bar;
                
                r_bar = 0.0;
                for k = 1:obj.n_nodes
                    f_k = U((1:4) + (k-1)*4);
                    q_k = f_k(1:4);
                    
                    p_k = -2.0 * (q_k(1) * q_bar(2:4) - q_bar(1) * q_k(2:4) - skew(q_k(2:4)) * q_bar(2:4));
                    Ua(ind_aux:ind_aux+2) = p_k;
                    ind_aux = ind_aux + 3;
                    
%                     r_k = sqrt(z2 - 0.25 * (p_k' * p_k));
                    r_k = q_bar(1) * q_k(1) - q_bar(2:4)' * q_k(2:4);
                    Ua(ind_aux) = r_k;
                    ind_aux = ind_aux + 1;
                    
                    r_bar = r_bar + shpfcn(k) * r_k;
                    
                    obj.F0(:, i) = obj.F0(:, i) + der_shpfcn(k) * p_k;
                end
                obj.F0(:, i) = obj.F0(:, i) / r_bar;
                ind_aux = ind_aux + 3;
                
                obj.J0(i) = 0.5;
                obj.F0(:, i) = obj.F0(:, i) / obj.J0(i);
            end
            
            loc_aux = loc_aux + ind_aux - 1;
        end
        
        function [R, Ra] = assemble(obj, U, Ua)
            
            R = zeros(3 * obj.n_nodes,1);
            Ra = zeros(obj.n_aux,1);
            
            ind_aux = 1;
            for i = 1:obj.n_gp
                eta = obj.loc_gp{obj.n_gp}(i);
                
                shpfcn = obj.shape_function{obj.n_nodes}(eta);
                der_shpfcn = obj.der_shape_function{obj.n_nodes}(eta) / obj.J0(i);
                
                q_bar = zeros(4,1);
                for k = 1:obj.n_nodes
                    f_k = U((1:4) + (k-1)*4);
                    q_bar = q_bar + shpfcn(k) * f_k(1:4);
                end
                z2 = q_bar' * q_bar;        

                r_bar = 0.0;
                p_k = zeros(3, obj.n_nodes);
                for k = 1:obj.n_nodes
                    f_k = U((1:4) + (k-1)*4);
                    q_k = f_k(1:4);
                    
                    p_k(:, k) = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = p_k(:, k) + 2.0 * (q_k(1) * q_bar(2:4) - q_bar(1) * q_k(2:4) - skew(q_k(2:4)) * q_bar(2:4));
                    ind_aux = ind_aux + 3;
                    
                    r_k = Ua(ind_aux);
                    Ra(ind_aux) = r_k * r_k - (z2 - 0.25 * (p_k(:, k)' * p_k(:, k)));
%                     Ra(ind_aux) = r_k * r_k - (q_bar(1) * q_k(1) - q_bar(2:4)' * q_k(2:4))^2;
                    ind_aux = ind_aux + 1;
                    
                    r_bar = r_bar + shpfcn(k) * r_k;
                end
                
                F = Ua(ind_aux:ind_aux+2);
                Ra(ind_aux:ind_aux+2) = r_bar * F;
                for k = 1:obj.n_nodes
                    Ra(ind_aux:ind_aux+2) = Ra(ind_aux:ind_aux+2) - der_shpfcn(k) * p_k(:, k);
                end
                ind_aux = ind_aux + 3;
                
                Q = zeros(3, 3 * obj.n_nodes);
                for k=1:obj.n_nodes
                    Q(:, (1:3) + (k-1)*3) = der_shpfcn(k) * eye(3) + shpfcn(k) * skew(F);
                end
                
                defo = F - obj.F0(:, i);
                stress = obj.beam_props.K * defo;
                
                R = R + Q' * (obj.weight_gp{obj.n_gp}(i) * obj.J0(i) * stress);
            end
            
        end
    end
end