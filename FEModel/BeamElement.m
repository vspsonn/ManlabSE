classdef BeamElement < Element

    properties
        beam_props

        n_gp = 2
        n_aux_gp

        J0
        G0
        R0
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
        function obj = BeamElement(list_node_numbers, beam_props, R0)
            obj = obj@Element();
            obj.number_of_nodes = size(list_node_numbers, 2);

            obj.loc_dof = zeros(1, 6 * obj.number_of_nodes);
            obj.loc_config = zeros(1, 7 * obj.number_of_nodes);
            for i = 1:obj.number_of_nodes
                obj.loc_dof((1:6) + (i-1) * 6) = [5:7 2:4] + 7 * (list_node_numbers(i)-1);
                obj.loc_config((1:7) + (i-1) * 7) = (1:7) + 7 * (list_node_numbers(i)-1);
            end

            obj.beam_props = beam_props;
            obj.R0 = R0;
        end

        function [rhs_aux, G, x_eta] = get_eq_gp(obj, F, Ua, i_gp, for_Ua0)
            n_cols = size(F, 2);

            rhs_aux = zeros(obj.n_aux_gp, n_cols);
            ind_aux = 1;

            eta = obj.loc_gp{obj.n_gp}(i_gp);
            shpfcn = obj.shape_function{obj.number_of_nodes}(eta);
            der_shpfcn = obj.der_shape_function{obj.number_of_nodes}(eta) / obj.J0(i_gp);
            G = zeros(6, n_cols);

            q_bar = zeros(4, n_cols);
            for k = 1:obj.number_of_nodes
                q_bar = q_bar + shpfcn(k) * F((1:4) + (k-1)*7, :);
            end

            a_R = q_bar(1, :) .* q_bar(1, :);
            a = Ua(ind_aux, :);
            rhs_aux(ind_aux, :) = a - a_R;
            if for_Ua0
                a = a_R;
            end
            ind_aux = ind_aux + 1;

            b = Ua(ind_aux:ind_aux+2, :);
            for i = 1:3
                b_R = q_bar(1, :) .* q_bar(i+1, :);
                rhs_aux(ind_aux, :) = b(i, :) - b_R;
                if for_Ua0
                    b(i, :) = b_R;
                end
                ind_aux = ind_aux + 1;
            end

            c = Ua(ind_aux:ind_aux+5, :);
            c_inds = [[1,1]; [2,2]; [3,3]; [1,2]; [1,3]; [2,3]];
            for i = 1:6
                ij = c_inds(i, :);
                c_R = q_bar(ij(1)+1, :) .* q_bar(ij(2)+1, :);
                rhs_aux(ind_aux, :) = c(i, :) - c_R;
                if for_Ua0
                    c(i, :) = c_R;
                end
                ind_aux = ind_aux + 1;
            end
            z2 = a + c(1, :) + c(2, :) + c(3, :);

            r_bar = zeros(1, n_cols);
            p_k = zeros(6, n_cols, obj.number_of_nodes);
            r_k = zeros(obj.number_of_nodes, n_cols);
            for k = 1:obj.number_of_nodes
                q_k = F((1:4) + (k-1)*7, :);

                p_k_R = - 2.0 * (q_k(1, :) .* q_bar(2:4, :) - q_bar(1, :) .* q_k(2:4, :) + SO3_cross(q_bar(2:4, :), q_k(2:4, :)));
                p_k(4:6, :, k) = Ua(ind_aux:ind_aux+2, :);
                rhs_aux(ind_aux:ind_aux+2, :) = p_k(4:6, :, k) - p_k_R;
                if for_Ua0
                    p_k(4:6, k) = p_k_R;
                end
                ind_aux = ind_aux + 3;

                r_k_R = z2 - 0.25 * (p_k(4, :, k) .* p_k(4, :, k) + p_k(5, :, k) .* p_k(5, :, k) + p_k(6, :, k) .* p_k(6, :, k));
                if for_Ua0
                    r_k(k) = sqrt(r_k_R);
                    rhs_aux(ind_aux) = - r_k(k);
                else
                    r_k(k, :) = Ua(ind_aux, :);
                    rhs_aux(ind_aux, :) = r_k(k, :) .* r_k(k, :) - r_k_R;
                end
                ind_aux = ind_aux + 1;

                r_bar = r_bar + shpfcn(k) * r_k(k, :);
            end

            if for_Ua0
                ind_aux_G46 = ind_aux;
                for k = 1:obj.number_of_nodes
                    G(4:6) = G(4:6) + der_shpfcn(k) * p_k(4:6, k);
                end
                G(4:6) = G(4:6) / r_bar;
            else
                G(4:6, :) = Ua(ind_aux:ind_aux+2, :);
                rhs_aux(ind_aux:ind_aux+2, :) = r_bar .* G(4:6, :);
                for k = 1:obj.number_of_nodes
                    rhs_aux(ind_aux:ind_aux+2, :) = rhs_aux(ind_aux:ind_aux+2, :) - der_shpfcn(k) * p_k(4:6, :, k);
                end
            end
            ind_aux = ind_aux + 3;

            x_eta_R = zeros(3, n_cols);
            Rp_k_R = zeros(3, n_cols);
            for k = 1:obj.number_of_nodes
                Rp_k_R(1, :) = 2 * ((-c(2, :) - c(3, :)) .* p_k(4, :, k) + (c(4, :) - b(3, :)) .* p_k(5, :, k) + (c(5, :) + b(2, :)) .* p_k(6, :, k));
                Rp_k_R(2, :) = 2 * ((c(4, :) + b(3, :)) .* p_k(4, :, k) + (-c(1, :) - c(3, :)) .* p_k(5, :, k) + (c(6, :) - b(1, :)) .* p_k(6, :, k));
                Rp_k_R(3, :) =  2 * ((c(5, :) - b(2, :)) .* p_k(4, :, k) + (c(6, :) + b(1, :)) .* p_k(5, :, k) + (-c(1, :) - c(2, :)) .* p_k(6, :, k));
                if for_Ua0
                    Rp_k = p_k(4:6, :, k) + Rp_k_R / z2;
                    rhs_aux(ind_aux:ind_aux+2) = - Rp_k;
                else
                    Rp_k = Ua(ind_aux:ind_aux+2, :);
                    rhs_aux(ind_aux:ind_aux+2, :) = z2 .* Rp_k - z2 .* p_k(4:6, :, k) - Rp_k_R;
                end
                ind_aux = ind_aux + 3;

                x_k = F((5:7) + (k-1)*7, :);
                x_eta_R = x_eta_R + shpfcn(k) * (r_k(k, :) .* x_k - 0.5 * SO3_cross(Rp_k, x_k));
            end
            if for_Ua0
                x_eta = x_eta_R / r_bar;
                rhs_aux(ind_aux:ind_aux+2) = - x_eta;
            else
                x_eta = Ua(ind_aux:ind_aux+2, :);
                rhs_aux(ind_aux:ind_aux+2, :) = r_bar .* x_eta - x_eta_R;
            end
            ind_aux = ind_aux + 3;

            b_bar_R = zeros(1, n_cols);
            d_k_R = zeros(3, n_cols);
            for k = 1:obj.number_of_nodes
                x_k = F((5:7) + (k-1)*7, :);
                Dx = x_k - x_eta;

                d_k_R(1, :) = 2 * ((-c(2, :) - c(3, :)) .* Dx(1, :) + (c(4, :) + b(3, :)) .* Dx(2, :) + (c(5, :) - b(2, :)) .* Dx(3, :));
                d_k_R(2, :) = 2 * ((c(4, :) - b(3, :)) .* Dx(1, :) + (-c(1, :) - c(3, :)) .* Dx(2, :) + (c(6, :) + b(1, :)) .* Dx(3, :));
                d_k_R(3, :) = 2 * ((c(5, :) + b(2, :)) .* Dx(1, :) + (c(6, :) - b(1, :)) .* Dx(2, :) + (-c(1, :) - c(2, :)) .* Dx(3, :));
                if for_Ua0
                    d_k = Dx + d_k_R / z2;
                    rhs_aux(ind_aux:ind_aux+2) = - d_k;
                else
                    d_k = Ua(ind_aux:ind_aux+2, :);
                    rhs_aux(ind_aux:ind_aux+2, :) = z2 .* d_k - z2 .* Dx - d_k_R;
                end
                ind_aux = ind_aux + 3;

                p_k_R = r_k(k, :) .* d_k - 0.5 * SO3_cross(p_k(4:6, :, k), d_k);
                p_k(1:3, :, k) = Ua(ind_aux:ind_aux+2, :);
                rhs_aux(ind_aux:ind_aux+2, :) = p_k(1:3, :, k) - p_k_R;
                if for_Ua0
                    p_k(1:3, :, k) = p_k_R;
                end
                ind_aux = ind_aux + 3;

                b_k = -0.25 * (d_k(1, :) .* p_k(4, :, k) + d_k(2, :) .* p_k(5, :, k) + d_k(3, :) .* p_k(6, :, k));
                b_bar_R = b_bar_R + shpfcn(k) * b_k;
            end

            if for_Ua0
                rhs_aux(ind_aux) = - b_bar_R;
                b_bar = b_bar_R;
            else
                b_bar = Ua(ind_aux, :);
                rhs_aux(ind_aux, :) = b_bar - b_bar_R;
            end
            ind_aux = ind_aux + 1;

            if for_Ua0
                for k = 1:obj.number_of_nodes
                    G(1:3) = G(1:3) + der_shpfcn(k) * p_k(1:3, :, k);
                end
                G(1:3) = G(1:3) - b_bar * G(4:6);
                G(1:3) = G(1:3) / r_bar;
                J0_ = norm(G(1:3));
                rhs_aux(ind_aux_G46:ind_aux_G46+2) = -G(4:6) / J0_;
                rhs_aux(ind_aux:ind_aux+2) = - G(1:3) / J0_;
            else
                G(1:3, :) = Ua(ind_aux:ind_aux+2, :);
                rhs_aux(ind_aux:ind_aux+2, :) = r_bar .* G(1:3, :) + b_bar .* G(4:6, :);
                for k = 1:obj.number_of_nodes
                    rhs_aux(ind_aux:ind_aux+2, :) = rhs_aux(ind_aux:ind_aux+2, :) - der_shpfcn(k) * p_k(1:3, :, k);
                end
            end

        end

        function QT = get_QT(obj, G, i_gp)
            eta = obj.loc_gp{obj.n_gp}(i_gp);
            shpfcn = obj.shape_function{obj.number_of_nodes}(eta);
            der_shpfcn = obj.der_shape_function{obj.number_of_nodes}(eta) / obj.J0(i_gp);

            n_cols = size(G, 2);
            QT = zeros(6 * obj.number_of_nodes, 6, n_cols);
            for j = 1:n_cols
                sG_theta_T = skew(G(4:6, j))';
                sG_u_T = skew(G(1:3, j))';
                for k = 1:obj.number_of_nodes
                    QT((1:3) + (k-1)*6, 1:3, j) = der_shpfcn(k) * eye(3) + shpfcn(k) * sG_theta_T;
                    QT((4:6) + (k-1)*6, 1:3, j) = shpfcn(k) * sG_u_T;
                    QT((4:6) + (k-1)*6, 4:6, j) = QT((1:3) + (k-1)*6, 1:3, j);
                end
            end
        end

        function [obj, Ua, loc_aux] = initialize(obj, ~, U, loc_aux)
            obj.n_aux_gp = 13 * obj.number_of_nodes + 20;
            obj.number_of_aux = obj.n_gp * obj.n_aux_gp;
            obj.loc_aux = loc_aux;

            F = U;
            for k = 1:obj.number_of_nodes
                if obj.R0{k}(1) == 1.0
                    continue;
                end
                F((1:4) + (k-1)*7) = SO3_composition(U((1:4) + (k-1)*7), obj.R0{k});
            end

            Ua = zeros(obj.number_of_aux, 1);
            obj.J0 = ones(1, obj.n_gp);
            obj.G0 = zeros(6, obj.n_gp);
            for i = 1:obj.n_gp
                [rhs_aux, G, ~] = get_eq_gp(obj, F, zeros(obj.n_aux_gp,1), i, true);

                obj.J0(i) = norm(G(1:3));
                obj.G0(:, i) = G / obj.J0(i);

                Ua((1:obj.n_aux_gp) + (i-1) * obj.n_aux_gp) = - rhs_aux;
            end

            loc_aux = loc_aux + obj.number_of_aux;
        end

        function [R, Ra] = assemble(obj, U, Ua)
            n_cols = size(U, 2);
            F = U;
            for k = 1:obj.number_of_nodes
                if obj.R0{k}(1) == 1.0
                    continue;
                end
                F((1:4) + (k-1)*7, :) = SO3_composition(U((1:4) + (k-1)*7, :), obj.R0{k});
            end

            R = zeros(6 * obj.number_of_nodes, n_cols);
            Ra = zeros(obj.number_of_aux, n_cols);
            for i = 1:obj.n_gp
                Ua_i = Ua((1:obj.n_aux_gp) + (i-1) * obj.n_aux_gp, :);
                [rhs_aux, G, ~] = get_eq_gp(obj, F, Ua_i, i, false);

                Ra((1:obj.n_aux_gp) + (i-1) * obj.n_aux_gp, :) = rhs_aux;

                deformation = G  - obj.G0(:, i);
                stress = obj.beam_props.K * deformation;

                stress = obj.weight_gp{obj.n_gp}(i) * obj.J0(i) * stress;
                QT =  obj.get_QT(G, i);
                if n_cols > 1
                    R = R + squeeze(sum(QT .* reshape(stress, 1, 6, n_cols), 2));
                else
                    R = R + QT * stress;
                end
            end

            for k = 1:obj.number_of_nodes
                if obj.R0{k}(1) == 1.0
                    continue;
                end
                R((1:3) + (k-1)*6, :) = SO3_RotateVector(obj.R0{k}, R((1:3) + (k-1)*6, :));
                R((4:6) + (k-1)*6, :) = SO3_RotateVector(obj.R0{k}, R((4:6) + (k-1)*6, :));
            end
        end

        function plot(obj, ~, U)
            if obj.number_of_nodes > 2
                X = zeros(3, obj.number_of_nodes);
                for k = 1:obj.number_of_nodes
                    X(:, k) = U((5:7) + (k-1)*7);
                end
            else
                F = zeros(7 * obj.number_of_nodes, 1);
                for k = 1:obj.number_of_nodes
                    F((1:4) + (k-1)*7) = SO3_composition(U((1:4) + (k-1)*7), obj.R0{k});
                    F((5:7) + (k-1)*7) = U((5:7) + (k-1)*7);
                end

                X = zeros(3, obj.n_gp+2);
                X(:,1) = U(5:7);
                for i = 1:obj.n_gp
                    [~, ~, x_eta] = get_eq_gp(obj, F, zeros(obj.n_aux_gp,1), i, true);
                    X(:, i+1) = x_eta;
                end
                X(:,end) = U(end-2:end);
            end

            %             for i = 1:size(X,2)-1
            %                 plot3([X(1,i) X(1,i+1)], [X(2,i) X(2,i+1)], [X(3,i) X(3,i+1)])
            %             end
            fnplt(cscvn(X))
        end
    end
end
