classdef ExternalForceElement < Element

    properties
        force_props
        lambda_dof
    end

    methods
        function obj = ExternalForceElement(node_number, force_props)
            obj.force_props = force_props;

            obj.loc_dof = zeros(1, 6);
            obj.loc_config = zeros(1, 7);
            obj.loc_dof = (2:7) + 7 * (node_number-1);
            obj.loc_config = (1:7) + 7 * (node_number-1);
        end

        function [obj, Ua, loc_aux] = initialize(obj, fe_model, U, loc_aux)
            obj.loc_config = [obj.loc_config fe_model.n_dof+1];


            if obj.force_props.follower
                obj.n_aux = 0;
                obj.loc_aux = [];
                Ua = [];
            else
                obj.n_aux = 6;
                obj.loc_aux = loc_aux;
                loc_aux = loc_aux + 6;
                Ua = zeros(6,1);
            end
        end

        function [R, Ra] = assemble(obj, U, Ua)
            F = obj.force_props.dir * U(8);

            if obj.force_props.follower
                Ra = zeros(0,1);

                R = zeros(6,1);
                R(1:3) = F(1:3);
                R(4:6) = F(4:6);
            else
                q0 = U(1);
                q = U(2:4);
                skq = skew(q);

                Ra = zeros(6,1);
                yu = Ua(1:3);
                Ra(1:3) = yu - 2.0 * skq * F(1:3);
                yr = Ua(4:6);
                Ra(4:6) = yr - 2.0 * skq * F(4:6);

                R = zeros(6,1);
                R(1:3) = F(1:3) - q0 * yu + skq * yu;
                R(4:6) = F(4:6) - q0 * yr + skq * yr;
            end
        end

    end

end