% global U Section Diagram   % Global variables to export point from the diagram.

%% parameters of the system

fe_model = FEModel(2);

beam_props.K = diag([1e8 1e8 1e6]);

n_beam_elem = 10;
beam_def = GeometryBeam2DDef(@deep_circular_arch_geometry, beam_props, n_beam_elem);
beam_def.number_of_nodes_per_beam = 3;
fe_model = fe_model.add_component(beam_def);

bc_node = ComponentPartDef(beam_def, 1);
fe_model = fe_model.add_bc_trans(bc_node);

bc2_node = ComponentPartDef(beam_def, -1);
fe_model = fe_model.add_bc_trans(bc2_node);
fe_model = fe_model.add_bc_rot(bc2_node);

ef_props.dir = [0; -1; 0];
ef_props.follower = false;
ef_node = ComponentPartDef(beam_def, 0.5);
ef_def = ExternalForceDef(ef_node, ef_props);
fe_model = fe_model.add_component(ef_def);

fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 1);
tracked_dof = tracked_dofs(4);

%% initialization of the system
sys=SystAQ(fe_model.n_dof_free,fe_model.n_aux,@fem_equations,@fem_point_display,@deep_circular_arch_global_display,fe_model,'vectorial');

%% Launch Manlab
Manlab('sys'          ,sys, ...        % description of your system
    'U0value'         ,fe_model.get_U0(0.0), ...         % starting point
    'displayvariables',[fe_model.n_dof_free+1 tracked_dof], ...    % MANLAB run
    'order'           ,6, ...     % order of the series
    'ANMthreshold'    ,1e-8, ...   % threshold for the domain of validity of the series
    'Amax_max'        ,1e6, ...    % maximum value of the domain of validity of the series
    'NRthreshold'     ,1e-7, ...   % threshold for Newton-Raphson (NR) corrections
    'NRitemax'        ,10, ...     % Maximum number of iteration of NR algorithm
    'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
    'NRmethod'        ,0, ...      % NR corrections on/[off]
    'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
    'PointDisplay'    ,0, ...      % Point display [on]/off
    'GlobalDisplay'   ,1, ...      % Global display [on]/off
    'StabilityCheck'  ,0, ...      % Stability computation on/[off]
    'StabTol'         ,1e-6);   % Stability tolerance
