global U Section Diagram  % Global variables to export point from the diagram.

%% parameters of the system

fe_model = FEModel(3);

beam_props.K = diag([10.0e6 4.167e6 4.167e6 83.34e4 83.34e4 83.34e4]);

n_beam_elem = 10;
beam_def = GeometryBeamDef(@bend45_geometry, beam_props, n_beam_elem);
beam_def.number_of_nodes_per_beam = 2;
fe_model = fe_model.add_component(beam_def);

bc_node = ComponentPartDef(beam_def, 1);
fe_model = fe_model.add_bc_trans(bc_node);
fe_model = fe_model.add_bc_rot(bc_node);

ef_props.dir = [0; 0; 1; 0; 0; 0];
ef_props.follower = false;
ef_node = ComponentPartDef(beam_def, -1);
ef_def = ExternalForceDef(ef_node, ef_props);
fe_model = fe_model.add_component(ef_def);

fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 1);
tracked_dof = tracked_dofs(7);

%% initialization of the system
sys=SystAQ(fe_model.n_dof_free,fe_model.n_aux,@fem_equations,@fem_point_display,@bend45_global_display,fe_model,'vectorial');

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
    'GlobalDisplay'   ,0, ...      % Global display [on]/off
    'StabilityCheck'  ,0, ...      % Stability computation on/[off]
    'StabTol'         ,1e-6);   % Stability tolerance
