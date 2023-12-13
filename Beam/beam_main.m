global U Section Diagram   % Global variables to export point from the diagram.

%% parameters of the system

fe_model = FEModel();

phi = 0.25 * pi;
f_root = [cos(phi/2); sin(phi/2)*[0.0; 0.0; 1.0]; [0.0; 0.0; 0.0]];
f_tip =  [cos(phi/2); sin(phi/2)*[0.0; 0.0; 1.0]; [1.0; 0.0; 0.0]];

beam_props.K = eye(6);

n_beam_elem = 3;
beam_def = StraightBeamDef(f_root, f_tip, beam_props, n_beam_elem);
beam_def.n_nodes_per_beam = 2;
fe_model = fe_model.add_component(beam_def);

bc_node = 1;
bc_def = BoundaryConditionDef(bc_node, 6);
fe_model = fe_model.add_component(bc_def);

ef_props.dir = [0; 0; 1; 0; 0; 0];
ef_props.follower = false;
ef_node = 4;
ef_def = ExternalForceDef(ef_node, ef_props);
fe_model = fe_model.add_component(ef_def);


fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

%% initialization of the system - Should not be changed
sys=SystAQ(fe_model.n_dof,fe_model.n_aux,@fem_equations,@point_display,@global_display,fe_model,'standard');

%% Launch Manlab
Manlab('sys'          ,sys, ...        % description of your system
    'U0value'         ,[fe_model.U0; 0; fe_model.Ua0], ...         % starting point
    'displayvariables',[fe_model.n_dof+1 fe_model.n_dof-6], ...    % MANLAB run
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
