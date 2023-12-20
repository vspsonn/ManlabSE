classdef StraightBeamDef < ComponentDef
    properties
        number_of_nodes_per_beam = 2
        
        frame_root
        frame_tip
        
        n_root = -1;
        n_tip = -1;
        
        beam_props
    end
    
    methods
        function obj = StraightBeamDef(frame_root, frame_tip, beam_props, number_of_elements)
            obj = obj@ComponentDef(3);
            
            obj.frame_root = frame_root;
            obj.frame_tip = frame_tip;
            obj.beam_props = beam_props;
            
            obj.number_of_elements = number_of_elements;
        end
        
        function [obj, U, FEModel_index_nodes] = mesh(obj, FEModel_index_nodes)
            obj.list_elements = cell(1, obj.number_of_elements);
            
            list_R0 = cell(obj.number_of_nodes_per_beam, 1);
            list_node_elem_numbers = zeros(1, obj.number_of_nodes_per_beam);
            
            number_of_nodes = obj.number_of_nodes_per_beam + (obj.number_of_elements-1) * (obj.number_of_nodes_per_beam-1);
            
            d = SE3_Log(SE3_composition(SE3_inverse(obj.frame_root), obj.frame_tip));
            
            number_of_nodes_to_mesh = number_of_nodes;
            if isnumeric(obj.n_root) && obj.n_root == -1
                FEModel_index_nodes = FEModel_index_nodes + 1;
                index_first_node = FEModel_index_nodes;
            else
                number_of_nodes_to_mesh = number_of_nodes_to_mesh - 1;
                index_first_node = obj.n_root.get_node_number();
            end
            if isnumeric(obj.n_tip) && obj.n_tip ~= -1
                number_of_nodes_to_mesh = number_of_nodes_to_mesh - 1;
            end
            
            U = zeros(7 * number_of_nodes_to_mesh, 1);
            for i = 1:number_of_nodes_to_mesh
                U(1 + (i-1) * 7) = 1;
            end
            
            list_R0{1} = obj.frame_root(1:4);
            list_node_elem_numbers(1) = index_first_node;
            obj.list_node_numbers = [index_first_node];
            if isnumeric(obj.n_root) && obj.n_root == -1
                index_node = 1;
                U(5:7) = obj.frame_root(5:7);
            else
                index_node = 0;
            end
            
            for i = 1:obj.number_of_elements
                for j = 1:(obj.number_of_nodes_per_beam-2)
                    alpha = ((i-1) + j / (obj.number_of_nodes_per_beam-1)) / obj.number_of_elements;
                    f_j = SE3_composition(obj.frame_root, SE3_Exp(alpha * d));
                    U((5:7) + index_node*7) = f_j(5:7);
                    list_R0{j+1} = f_j(1:4);
                    
                    index_node = index_node + 1;
                    FEModel_index_nodes = FEModel_index_nodes + 1;
                    list_node_elem_numbers(j+1) = FEModel_index_nodes;
                end
                if i == obj.number_of_elements
                    list_R0{end} = obj.frame_tip(1:4);
                    if isnumeric(obj.n_tip) && obj.n_tip == -1
                        U(end-2:end) = obj.frame_tip(5:7);
                        
                        FEModel_index_nodes = FEModel_index_nodes + 1;
                        list_node_elem_numbers(obj.number_of_nodes_per_beam) = FEModel_index_nodes;
                    else
                        list_node_elem_numbers(obj.number_of_nodes_per_beam) = obj.n_tip.get_node_number();
                    end
                else
                    alpha = i / obj.number_of_elements;
                    f_i = SE3_composition(obj.frame_root, SE3_Exp(alpha * d));
                    U((5:7) + index_node*7) = f_i(5:7);
                    list_R0{end} = f_i(1:4);
                    
                    index_node = index_node + 1;
                    FEModel_index_nodes = FEModel_index_nodes + 1;
                    list_node_elem_numbers(obj.number_of_nodes_per_beam) = FEModel_index_nodes;
                end
                
                obj.list_elements{i} = BeamElement(list_node_elem_numbers, obj.beam_props, list_R0);
                obj.list_elements{i}.n_gp = obj.number_of_nodes_per_beam-1;
                
                obj.list_node_numbers = [obj.list_node_numbers list_node_elem_numbers(2:end)];
                list_node_elem_numbers(1) = list_node_elem_numbers(end);
                list_R0{1} = list_R0{end};
            end
        end
        
    end
end