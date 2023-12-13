classdef BeamElement < Element
    
    properties
        beam_props
        
        n_gp = 2
        J0
        G0
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
        function obj = BeamElement(list_node_numbers, beam_props)
            obj.n_nodes = size(list_node_numbers, 2);
            
            obj.loc_dof = zeros(1, 6 * obj.n_nodes);
            obj.loc_config = zeros(1, 7 * obj.n_nodes);
            for i = 1:obj.n_nodes
                obj.loc_dof((1:6) + (i-1) * 6) = (2:7) + 7 * (list_node_numbers(i)-1);
                obj.loc_config((1:7) + (i-1) * 7) = (1:7) + 7 * (list_node_numbers(i)-1);
            end
            
            obj.beam_props = beam_props;
        end
        
        function [obj, Ua, loc_aux] = initialize(obj, fe_model, U, loc_aux)
            
            obj.n_aux = obj.n_gp * (19 * obj.n_nodes + 11);
            obj.loc_aux = loc_aux;
            
            Ua = zeros(obj.n_aux, 1);
            obj.J0 = zeros(1, obj.n_gp);
            obj.G0 = zeros(6, obj.n_gp);
            
            ind_aux = 1;
            for i = 1:obj.n_gp
                eta = obj.loc_gp{obj.n_gp}(i);
                
                shpfcn = obj.shape_function{obj.n_nodes}(eta);
                der_shpfcn = obj.der_shape_function{obj.n_nodes}(eta);
                
                q_bar = zeros(4,1);
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    q_bar = q_bar + shpfcn(k) * F_k(1:4);
                end
                z = sqrt(q_bar'*q_bar);
                Ua(ind_aux) = z;
                ind_aux = ind_aux + 1;
                
                r_bar = 0.0;
                p_k = zeros(6, obj.n_nodes);
                r_k = zeros(1, obj.n_nodes);
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    q_k = F_k(1:4);
                    
                    p_k(4:6, k) = -2.0 / z * (q_k(1) * q_bar(2:4) - q_bar(1) * q_k(2:4) - skew(q_k(2:4)) * q_bar(2:4));
                    Ua(ind_aux:ind_aux+2) = p_k(4:6, k);
                    ind_aux = ind_aux + 3;
                    
                    r_k(k) = sqrt(1. - 0.25 * (p_k(4:6, k)'*p_k(4:6, k)));
                    Ua(ind_aux) = r_k(k);
                    ind_aux = ind_aux + 1;
                    
                    r_bar = r_bar + shpfcn(k) * r_k(k);
                    
                    obj.G0(4:6, i) = obj.G0(4:6, i) + der_shpfcn(k) * p_k(4:6, k);
                end
                ind_aux_G0_rot = ind_aux;
                obj.G0(4:6, i) = obj.G0(4:6, i) / r_bar;
                ind_aux = ind_aux + 3;
                
                x_eta = zeros(3,1);
                for k = 1:obj.n_nodes
                    yk = skew(q_bar(2:4)) * p_k(4:6, k) / z;
                    Ua(ind_aux:ind_aux+2) = yk;
                    ind_aux = ind_aux + 3;
                    
                    Rp_k = p_k(4:6, k) + 2.0 / z * (q_bar(1) * yk + skew(q_bar(2:4)) * yk);
                    Ua(ind_aux:ind_aux+2) = Rp_k;
                    ind_aux = ind_aux + 3;
                    
                    F_k = U((1:7) + (k-1)*7);
                    x_k = F_k(5:7);
                    x_eta = x_eta + shpfcn(k) / r_bar * (r_k(k) * x_k - 0.5* skew(Rp_k) * x_k);
                end
                Ua(ind_aux:ind_aux+2) = x_eta;
                ind_aux = ind_aux + 3;
                
                b_bar = 0.0;
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    x_k = F_k(5:7);
                    Dx = x_k - x_eta;
                    
                    u_k = skew(q_bar(2:4)) * Dx / z;
                    Ua(ind_aux:ind_aux+2) = u_k;
                    ind_aux = ind_aux + 3;
                    
                    d_k = Dx - 2.0 / z * (q_bar(1) * u_k - skew(q_bar(2:4)) * u_k);
                    Ua(ind_aux:ind_aux+2) = d_k;
                    ind_aux = ind_aux + 3;
                    
                    p_k(1:3, k) = r_k(k) * d_k - 0.5 * skew(p_k(4:6, k)) * d_k;
                    Ua(ind_aux:ind_aux+2) = p_k(1:3, k);
                    ind_aux = ind_aux + 3;
                    
                    %b_k = -0.25 / r_k(k) * (p_k(1:3, k)' * p_k(4:6, k));
                    b_k = -0.25 * d_k' * p_k(4:6, k);
                    b_bar = b_bar + shpfcn(k) * b_k;
                    
                    obj.G0(1:3, i) = obj.G0(1:3, i) + der_shpfcn(k) * p_k(1:3, k);
                end
                Ua(ind_aux) = b_bar;
                ind_aux = ind_aux + 1;

                obj.G0(1:3, i) = obj.G0(1:3, i) - b_bar * obj.G0(4:6, i);
                obj.G0(1:3, i) = obj.G0(1:3, i) / r_bar;
                
                obj.J0(i) = norm(obj.G0(1:3, i));
                obj.G0(:, i) = obj.G0(:, i) / obj.J0(i);
                
                Ua(ind_aux:ind_aux+2) = obj.G0(1:3, i);
                ind_aux = ind_aux + 3;
                
                Ua(ind_aux_G0_rot:ind_aux_G0_rot+2) = obj.G0(4:6, i);
            end
            
            loc_aux = loc_aux + ind_aux - 1;
        end
        
        function [R, Ra] = assemble(obj, U, Ua)
            
            R = zeros(6*obj.n_nodes,1);
            Ra = zeros(obj.n_aux,1);
            
            ind_aux = 1;
            for i = 1:obj.n_gp
                eta = obj.loc_gp{obj.n_gp}(i);
                
                shpfcn = obj.shape_function{obj.n_nodes}(eta);
                der_shpfcn = obj.der_shape_function{obj.n_nodes}(eta) / obj.J0(i);
                
                q_bar = zeros(4,1);
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    q_bar = q_bar + shpfcn(k) * F_k(1:4);
                end
                
                z = Ua(ind_aux);
                Ra(ind_aux) = z * z - q_bar'*q_bar;
                ind_aux = ind_aux + 1;
                
                r_bar = 0.0;
                p_k = zeros(6, obj.n_nodes);
                r_k = zeros(1, obj.n_nodes);
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    q_k = F_k(1:4);
                    
                    p_k(4:6, k) = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = z * p_k(4:6, k) + 2.0 * (q_k(1) * q_bar(2:4) - q_bar(1) * q_k(2:4) - skew(q_k(2:4)) * q_bar(2:4));
                    ind_aux = ind_aux + 3;
                    
                    r_k(k) = Ua(ind_aux);
                    Ra(ind_aux) = r_k(k) * r_k(k) - (1. - 0.25 * (p_k(4:6, k)' * p_k(4:6, k)));
                    ind_aux = ind_aux + 1;
                    
                    r_bar = r_bar + shpfcn(k) * r_k(k);
                end
                
                G = zeros(6, 1);
                G(4:6) = Ua(ind_aux:ind_aux+2);
                Ra(ind_aux:ind_aux+2) = r_bar * G(4:6);
                for k=1:obj.n_nodes
                    Ra(ind_aux:ind_aux+2) = Ra(ind_aux:ind_aux+2) - der_shpfcn(k) * p_k(4:6, k);
                end
                ind_aux = ind_aux + 3;
                
                Rp_k = zeros(3, obj.n_nodes);
                for k = 1:obj.n_nodes
                    yk = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = z * yk - skew(q_bar(2:4)) * p_k(4:6, k);
                    ind_aux = ind_aux + 3;
                    
                    Rp_k(:, k) = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = z * Rp_k(:, k) - z * p_k(4:6, k) - 2.0 * (q_bar(1) * yk + skew(q_bar(2:4)) * yk);
                    ind_aux = ind_aux + 3;
                end
                
                x_eta = Ua(ind_aux:ind_aux+2);
                Ra(ind_aux:ind_aux+2) = r_bar * x_eta;
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    x_k = F_k(5:7);
                    Ra(ind_aux:ind_aux+2) = Ra(ind_aux:ind_aux+2) - shpfcn(k) * (r_k(k) * x_k - 0.5* skew(Rp_k(:, k)) * x_k);
                end
                ind_aux = ind_aux + 3;
                
                b_bar_comp = 0.0;
                for k = 1:obj.n_nodes
                    F_k = U((1:7) + (k-1)*7);
                    x_k = F_k(5:7);
                    Dx = x_k - x_eta;
                    
                    u_k = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = z * u_k - skew(q_bar(2:4)) * Dx;
                    ind_aux = ind_aux + 3;
                    
                    d_k = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = z * d_k - z * Dx + 2.0 * (q_bar(1) * u_k - skew(q_bar(2:4)) * u_k);
                    ind_aux = ind_aux + 3;
                    
                    p_k(1:3, k) = Ua(ind_aux:ind_aux+2);
                    Ra(ind_aux:ind_aux+2) = p_k(1:3, k) - (r_k(k) * d_k - 0.5 * skew(p_k(4:6, k)) * d_k);
                    ind_aux = ind_aux + 3;
                    
                    b_k = 0.25 * d_k' * p_k(4:6, k);
                    b_bar_comp = b_bar_comp + shpfcn(k) * b_k;
                end
                b_bar = Ua(ind_aux);
                Ra(ind_aux) = b_bar + b_bar_comp;
                ind_aux = ind_aux + 1;
                
                G(1:3) = Ua(ind_aux:ind_aux+2);
                Ra(ind_aux:ind_aux+2) = r_bar * G(1:3) + b_bar * G(4:6);
                for k=1:obj.n_nodes
                    Ra(ind_aux:ind_aux+2) = Ra(ind_aux:ind_aux+2) - der_shpfcn(k) * p_k(1:3, k);
                end
                ind_aux = ind_aux + 3;
                
                Q = zeros(6, 6*obj.n_nodes);
                for k=1:obj.n_nodes
                    Q(:, (1:6) + (k-1)*6) = der_shpfcn(k) * eye(6) + shpfcn(k) * [skew(G(4:6)) skew(G(1:3)); zeros(3) skew(G(4:6))];
                end
                
                defo = G - obj.G0(:, i);
                stress = obj.beam_props.K * defo;
                
                R = R + Q' * (obj.weight_gp{obj.n_gp}(i) * obj.J0(i) * stress);
            end
            
        end
        
        function plot(obj, fig, U)
            x_A = U(5:7);
            x_B = U(end-2:end);
            plot3([x_A(1) x_B(1)], [x_A(2) x_B(2)], [x_A(3) x_B(3)])            
        end
    end
end