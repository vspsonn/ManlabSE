global U Section Diagram   % Global variables to export point from the diagram.

%% parameters of the system

fe_model = FEModel(2);

E = 7.2e6;
nu = 0.3;
G = E / (2.0*(1.0 + nu));
A = 6.0;
I = 2.0;

beam_props.K = diag([E*A, G*A, E*I]);
beam_length = 120;

phi = 0.5 * pi;
q1 = [cos(phi/2); sin(phi/2)];
f1_root = [q1; [0.0; 0.0]];
f1_tip =  [q1; [0.0; beam_length]];

n_beam_elem = 10;
beam_def = StraightBeam2DDef(f1_root, f1_tip, beam_props, n_beam_elem);
beam_def.number_of_nodes_per_beam = 2;
fe_model = fe_model.add_component(beam_def);

f2_root = [[1; 0]; [0.0; beam_length]];
f2_tip =  [[1; 0]; [beam_length; beam_length]];

n_beam_elem = 10;
beam2_def = StraightBeam2DDef(f2_root, f2_tip, beam_props, n_beam_elem);
beam2_def.number_of_nodes_per_beam = 2;
beam2_def.n_root = ComponentPartDef(beam_def, -1);
fe_model = fe_model.add_component(beam2_def);

bc_node = ComponentPartDef(beam_def, 1);
fe_model = fe_model.add_bc_trans(bc_node);

bc2_node = ComponentPartDef(beam2_def, -1);
fe_model = fe_model.add_bc_trans(bc2_node);

ef_props.dir = [0; -1; 0];
ef_props.follower = false;
ef_node = ComponentPartDef(beam2_def, 0.2);
ef_def = ExternalForceDef(ef_node, ef_props);
fe_model = fe_model.add_component(ef_def);

ref_node = ComponentPartDef(beam2_def, 0);

fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

ref_dofs = fe_model.get_nodal_dofs(ref_node, 1);
tracked_dof = ref_dofs(3);

%% initialization of the system
sys=SystAQ(fe_model.n_dof_free,fe_model.n_aux,@fem_equations,@fem_point_display,@hinged_frame_global_display,fe_model,'vectorial');

%% Launch Manlab
Manlab('sys'          ,sys, ...        % description of your system
    'U0value'         ,fe_model.get_U0(0.0), ...         % starting point
    'displayvariables',[tracked_dof fe_model.n_dof_free+1], ...    % MANLAB run
    'order'           ,11, ...     % order of the series
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
