global U Section Diagram   % Global variables to export point from the diagram.

%% parameters of the system

fe_model = FEModel(2);

phi = 0. * pi;
q = [cos(phi/2); sin(phi/2)];
f_root = [q; [0.0; 0.0]];
f_tip =  [q; SO2_RotateVector(q, [1.0; 0.0])];

beam_props.K = eye(3);

n_beam_elem = 3;
beam_def = StraightBeam2DDef(f_root, f_tip, beam_props, n_beam_elem);
beam_def.number_of_nodes_per_beam = 2;
fe_model = fe_model.add_component(beam_def);

bc_node = ComponentPartDef(beam_def, 1);
fe_model = fe_model.add_bc_trans(bc_node);
fe_model = fe_model.add_bc_rot(bc_node);

ef_props.dir = [0; 0; 1];
ef_props.follower = true;
ef_node = ComponentPartDef(beam_def, -1);
ef_def = ExternalForceDef(ef_node, ef_props);
fe_model = fe_model.add_component(ef_def);

fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

tracked_dofs = fe_model.get_nodal_dofs(ef_node, 1);
tracked_dof = tracked_dofs(4);

%% initialization of the system - Should not be changed
sys=SystAQ(fe_model.n_dof_free,fe_model.n_aux,@fem_equations,@fem_point_display,@global_display,fe_model,'vectorial');

%% Launch Manlab
Manlab('sys'          ,sys, ...        % description of your system
    'U0value'         ,fe_model.get_U0(0.0), ...         % starting point
    'displayvariables',[fe_model.n_dof_free+1 tracked_dof], ...    % MANLAB run
    'order'           ,5, ...     % order of the series
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

% With 'displayvariables',[neq+1 1] ; the projected bifurcation diagram
% (figure 2) will represents the first variable with respect to the
% continuation parameter lambda.
%
%               The stability requires to write the system of equations in
%               the explicit ODE form X' = f(X).
%               The residue function is then R(X) = f(X) = 0.
%               In all other cases, it gives wrong results.
