classdef FEModel
    
    properties
        list_components
        
        n_dof
        U0
        
        n_aux
        Ua0
    end
    
    methods
        function obj = FEModel()
            obj.list_components = cell(1,0);
        end
        
        function obj = add_component(obj, component)
            obj.list_components{size(obj.list_components, 2)+1} = component;
        end
        
        function obj = mesh(obj)
            obj.U0 = [];
            ind_nodes = 1;
            for i = 1:size(obj.list_components, 2)
                [obj.list_components{i}, U0i, ind_nodes] = obj.list_components{i}.mesh(ind_nodes);
                obj.U0 = [obj.U0; U0i];
            end
            
            obj.n_dof = size(obj.U0, 1);
        end
        
        function obj = initialize(obj)
            obj.Ua0 = [];
            ind_aux = 1;
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elem, 2)
                    elem = obj.list_components{i}.list_elem{j};
                    [obj.list_components{i}.list_elem{j}, Ua0_j, ind_aux] = elem.initialize(obj, obj.U0(elem.loc_config), ind_aux);
                    obj.Ua0 = [obj.Ua0; Ua0_j];
                end
            end
            
            obj.n_aux = size(obj.Ua0, 1);
        end
        
        function [R, Ra] = assemble(obj, Uf)
            
            U = Uf(1:obj.n_dof+1);
            Ua = Uf(obj.n_dof+2:end);
            
            R = zeros(obj.n_dof, 1);
            Ra = zeros(obj.n_aux, 1);
            
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elem, 2)
                    elem = obj.list_components{i}.list_elem{j};
                    loc_aux = elem.loc_aux + (0:(elem.n_aux-1));
                    [R_i, Ra_i] = elem.assemble(U(elem.loc_config), Ua(loc_aux));
                    R(elem.loc_dof) = R(elem.loc_dof) + R_i;
                    Ra(loc_aux) = Ra_i;
                end
            end
            
            Moment = U(obj.n_dof+1);
            %R(obj.n_dof-9) = R(obj.n_dof-9) + Moment;
            R(7) = R(7) + Moment;
        end
        
        function plot(obj, fig, U)
            for i = 1:size(obj.list_components, 2)
                for j = 1:size(obj.list_components{i}.list_elem, 2)
                    elem = obj.list_components{i}.list_elem{j};
                    elem.plot(fig, U(elem.loc_config));
                end
            end
        end
        
    end
end