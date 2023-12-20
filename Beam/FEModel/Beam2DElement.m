classdef Beam2DElement < Element

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
        function obj = Beam2DElement(list_node_numbers, beam_props, R0)
            obj = obj@Element();
            obj.number_of_nodes = size(list_node_numbers, 2);

            obj.loc_dof = zeros(1, 3 * obj.number_of_nodes);
            obj.loc_config = zeros(1, 4 * obj.number_of_nodes);
            for i = 1:obj.number_of_nodes
                obj.loc_dof((1:3) + (i-1) * 3) = [3:4 2] + 4 * (list_node_numbers(i)-1);
                obj.loc_config((1:4) + (i-1) * 4) = (1:4) + 4 * (list_node_numbers(i)-1);
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
            G = zeros(3, n_cols);

            q_bar = zeros(2, n_cols);
            for k = 1:obj.number_of_nodes
                q_bar = q_bar + shpfcn(k) * F((1:2) + (k-1)*4, :);
            end

            a_R = q_bar(1, :) .* q_bar(1, :);
            a = Ua(ind_aux, :);
            rhs_aux(ind_aux, :) = a - a_R;
            if for_Ua0
                a = a_R;
            end
            ind_aux = ind_aux + 1;

            b_R = q_bar(1, :) .* q_bar(2, :);
            b = Ua(ind_aux, :);
            rhs_aux(ind_aux, :) = b - b_R;
            if for_Ua0
                b = b_R;
            end
            ind_aux = ind_aux + 1;

            c_R = q_bar(2, :) .* q_bar(2, :);
            c = Ua(ind_aux, :);
            rhs_aux(ind_aux, :) = c - c_R;
            if for_Ua0
                c = c_R;
            end
            ind_aux = ind_aux + 1;

            z2 = a + c;

            r_bar = zeros(1, n_cols);
            p_k = zeros(3, n_cols, obj.number_of_nodes);
            r_k = zeros(obj.number_of_nodes, n_cols);
            for k = 1:obj.number_of_nodes
                q_k = F((1:2) + (k-1)*4, :);

                p_k_R = - 2.0 * (q_k(1, :) .* q_bar(2, :) - q_bar(1, :) .* q_k(2, :));
                p_k(3, :, k) = Ua(ind_aux, :);
                rhs_aux(ind_aux, :) = p_k(3, :, k) - p_k_R;
                if for_Ua0
                    p_k(3, k) = p_k_R;
                end
                ind_aux = ind_aux + 1;

                r_k_R = z2 - 0.25 * p_k(3, :, k) .* p_k(3, :, k);
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
                ind_aux_G3 = ind_aux;
                for k = 1:obj.number_of_nodes
                    G(3) = G(3) + der_shpfcn(k) * p_k(3, k);
                end
                G(3) = G(3) / r_bar;
            else
                G(3, :) = Ua(ind_aux, :);
                rhs_aux(ind_aux, :) = r_bar .* G(3, :);
                for k = 1:obj.number_of_nodes
                    rhs_aux(ind_aux, :) = rhs_aux(ind_aux, :) - der_shpfcn(k) * p_k(3, :, k);
                end
            end
            ind_aux = ind_aux + 1;

            x_eta_R = zeros(2, n_cols);
            for k = 1:obj.number_of_nodes
                x_k = F((3:4) + (k-1)*4, :);
                x_eta_R = x_eta_R + shpfcn(k) * (r_k(k, :) .* x_k + (0.5*p_k(3, :, k)) .* [x_k(2, :); -x_k(1, :)]);
            end
            if for_Ua0
                x_eta = x_eta_R / r_bar;
                rhs_aux(ind_aux:ind_aux+1) = - x_eta;
            else
                x_eta = Ua(ind_aux:ind_aux+1, :);
                rhs_aux(ind_aux:ind_aux+1, :) = r_bar .* x_eta - x_eta_R;
            end
            ind_aux = ind_aux + 2;

            for k = 1:obj.number_of_nodes
                x_k = F((3:4) + (k-1)*4, :);
                Dx = x_k - x_eta;

                d_k_R = (-2.0*b) .* [-Dx(2, :); Dx(1, :)] - (2.0*c) .* Dx;
                if for_Ua0
                    d_k = Dx + d_k_R / z2;
                    rhs_aux(ind_aux:ind_aux+1) = - d_k;
                else
                    d_k = Ua(ind_aux:ind_aux+1, :);
                    rhs_aux(ind_aux:ind_aux+1, :) = z2 .* d_k - z2 .* Dx - d_k_R;
                end
                ind_aux = ind_aux + 2;

                p_k_R = r_k(k, :) .* d_k + (0.5*p_k(3, :, k)) .* [d_k(2, :); -d_k(1, :)];
                p_k(1:2, :, k) = Ua(ind_aux:ind_aux+1, :);
                rhs_aux(ind_aux:ind_aux+1, :) = p_k(1:2, :, k) - p_k_R;
                if for_Ua0
                    p_k(1:2, :, k) = p_k_R;
                end
                ind_aux = ind_aux + 2;
            end

            if for_Ua0
                for k = 1:obj.number_of_nodes
                    G(1:2) = G(1:2) + der_shpfcn(k) * p_k(1:2, :, k);
                end
                G(1:2) = G(1:2) / r_bar;
                J0_ = norm(G(1:2));
                rhs_aux(ind_aux_G3) = - G(3) / J0_;
                rhs_aux(ind_aux:ind_aux+1) = - G(1:2) / J0_;
            else
                G(1:2, :) = Ua(ind_aux:ind_aux+1, :);
                rhs_aux(ind_aux:ind_aux+1, :) = r_bar .* G(1:2, :);
                for k = 1:obj.number_of_nodes
                    rhs_aux(ind_aux:ind_aux+1, :) = rhs_aux(ind_aux:ind_aux+1, :) - der_shpfcn(k) * p_k(1:2, :, k);
                end
            end

        end

        function QT = get_QT(obj, G, i_gp)
            eta = obj.loc_gp{obj.n_gp}(i_gp);
            shpfcn = obj.shape_function{obj.number_of_nodes}(eta);
            der_shpfcn = obj.der_shape_function{obj.number_of_nodes}(eta) / obj.J0(i_gp);

            n_cols = size(G, 2);
            QT = zeros(3 * obj.number_of_nodes, 3, n_cols);
            for j = 1:n_cols
                G_hat_T = [0 -G(3, j) G(2, j); G(3, j) 0 -G(1, j)]';
                for k = 1:obj.number_of_nodes
                    QT((1:3) + (k-1)*3, 1:3, j) = der_shpfcn(k) * eye(3);
                    QT((1:3) + (k-1)*3, 1:2, j) = QT((1:3) + (k-1)*3, 1:2, j) + shpfcn(k) * G_hat_T;
                end
            end
        end

        function [obj, Ua, loc_aux] = initialize(obj, ~, U, loc_aux)
            obj.n_aux_gp = 6 * obj.number_of_nodes + 8;
            obj.number_of_aux = obj.n_gp * obj.n_aux_gp;
            obj.loc_aux = loc_aux;

            F = U;
            for k = 1:obj.number_of_nodes
                if obj.R0{k}(1) == 1.0
                    continue;
                end
                F((1:2) + (k-1)*4) = SO2_composition(U((1:2) + (k-1)*4), obj.R0{k});
            end

            Ua = zeros(obj.number_of_aux, 1);
            obj.J0 = ones(1, obj.n_gp);
            obj.G0 = zeros(3, obj.n_gp);
            for i = 1:obj.n_gp
                [rhs_aux, G, ~] = get_eq_gp(obj, F, zeros(obj.n_aux_gp,1), i, true);

                obj.J0(i) = norm(G(1:2));
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
                F((1:2) + (k-1)*4, :) = SO2_composition(U((1:2) + (k-1)*4, :), obj.R0{k});
            end

            R = zeros(3 * obj.number_of_nodes, n_cols);
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
                    R = R + squeeze(sum(QT .* reshape(stress, 1, 3, n_cols), 2));
                else
                    R = R + QT * stress;
                end
            end

            for k = 1:obj.number_of_nodes
                if obj.R0{k}(1) == 1.0
                    continue;
                end
                R((1:2) + (k-1)*3, :) = SO2_RotateVector(obj.R0{k}, R((1:2) + (k-1)*3, :));
            end
        end

        function plot(obj, ~, U)
            if obj.number_of_nodes > 2
                X = zeros(2, obj.number_of_nodes);
                for k = 1:obj.number_of_nodes
                    X(:, k) = U((3:4) + (k-1)*4);
                end
            else
                F = zeros(4*obj.number_of_nodes, 1);
                for k = 1:obj.number_of_nodes
                    F((1:2) + (k-1)*4) = SO2_composition(U((1:2) + (k-1)*4), obj.R0{k});
                    F((3:4) + (k-1)*4) = U((3:4) + (k-1)*4);
                end

                X = zeros(2, obj.n_gp+2);
                X(:,1) = U(3:4);
                for i = 1:obj.n_gp
                    [~, ~, x_eta] = get_eq_gp(obj, F, zeros(obj.n_aux_gp,1), i, true);
                    X(:, i+1) = x_eta;
                end
                X(:,end) = U(end-1:end);
            end

            %             for i = 1:size(X,2)-1
            %                 plot([X(1,i) X(1,i+1)], [X(2,i) X(2,i+1)])
            %             end
            fnplt(cscvn(X))
        end
    end
end
