classdef StraightBeamDef < ComponentDef
    properties
        n_nodes_per_beam = 2
        
        f_root
        f_tip
        
        n_root = -1
        n_tip = -1
        
        beam_props
    end
    
    methods
        function obj = StraightBeamDef(f_root, f_tip, beam_props, n_elem)
            obj = obj@ComponentDef();
            
            obj.f_root = f_root;
            obj.f_tip = f_tip;
            obj.beam_props = beam_props;
            
            obj.n_elem = n_elem;
        end
        
        function [obj, U, ind_nodes] = mesh(obj, ind_nodes)
            n_nodes = obj.n_nodes_per_beam + (obj.n_elem-1) * (obj.n_nodes_per_beam-1);
            obj.list_elem = cell(1, obj.n_elem + n_nodes);
            
            if obj.n_elem == 0
                return
            end
            
            d = SE3_Log(SE3_composition(SE3_inverse(obj.f_root), obj.f_tip));
            
            U = zeros(7 * n_nodes, 1);
            U(1:7) = obj.f_root;
            
            list_nodes = zeros(1, obj.n_nodes_per_beam);
            ind_nodes0 = ind_nodes - 1;
            list_nodes(1) = ind_nodes;
            for i = 1:obj.n_elem
                for j = 1:(obj.n_nodes_per_beam-2)
                    alpha = ((i-1) + j / (obj.n_nodes_per_beam-1)) / obj.n_elem;
                    U((1:7) + ind_nodes*7) = SE3_composition(obj.f_root, SE3_Exp(alpha * d));
                    
                    ind_nodes = ind_nodes + 1;
                    list_nodes(j+1) = ind_nodes;
                end
                if i == obj.n_elem
                    U((1:7) + ind_nodes*7) = obj.f_tip;
                else
                    alpha = i / obj.n_elem;
                    U((1:7) + ind_nodes*7) = SE3_composition(obj.f_root, SE3_Exp(alpha * d));
                end
                ind_nodes = ind_nodes + 1;
                list_nodes(obj.n_nodes_per_beam) = ind_nodes;
                
                obj.list_elem{i} = Beam2Element(list_nodes, obj.beam_props);
                obj.list_elem{i}.n_gp = 1;
                
                list_nodes(1) = list_nodes(obj.n_nodes_per_beam);
            end
            
            for i = 1:n_nodes
                obj.list_elem{obj.n_elem + i} = UnitQuatConstraint(ind_nodes0 + i, 6);
            end
            
        end
        
    end
end