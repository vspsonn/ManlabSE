global U Section Diagram   % Global variables to export point from the diagram.

%% parameters of the system

f_root = [1.0; zeros(3,1)];
f_tip = [1.0; zeros(3,1)];

beam_props.K = eye(3);

beam_def = StraightBeamFlexDef(f_root, f_tip, beam_props, 1);
beam_def.n_nodes_per_beam = 2;

bc_def = BoundaryConditionDef(1, 3);

fe_model = FEModel();
fe_model = fe_model.add_component(beam_def);
fe_model = fe_model.add_component(bc_def);

fe_model = fe_model.mesh();
fe_model = fe_model.initialize();

%% initialization of the system - Should not be changed
sys=SystAQ(fe_model.n_dof,fe_model.n_aux,@fem_equations,@point_display,@global_display,fe_model,'standard');

%% Launch Manlab
Manlab('sys'          ,sys, ...        % description of your system
    'U0value'         ,[fe_model.U0; 0; fe_model.Ua0], ...         % starting point
    'displayvariables',[fe_model.n_dof+1 fe_model.n_dof-4], ...    % MANLAB run
    'order'           ,20, ...     % order of the series
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
