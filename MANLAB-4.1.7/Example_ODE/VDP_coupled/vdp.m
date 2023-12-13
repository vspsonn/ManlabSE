% Launch Van der Pol problem
%
% Les equations sont : [i=1,2] pi''+ wi/Qi * pi' + wi*wi * pi = Fi*[A+B(p1+p2)+C(p1+p2)^2](p1'+p2')    -    Van der Pol couplés
%
% On reformule le système :
%
% p1' = q1
% q1' = - w1^2*p1+a1*q1+a2*lambda*(q1+q2)+a3*lambda*w+a4*lambda*v
% p2' = q2
% q2' =- w1^2*p2+b1*q2+b2*lambda*(q1+q2)+b3*lambda*w+b4*lambda*v
%  0   =r - (p1 + p2)^2
%  0   =w - (p1+p2)*(q1+q2)
%  0   =v- r*(q1+q2)
%
%    z = [ p1 ; q1 ; p2 ; q2 ; r ; w ; v ]  (nz_tot=7)  % unknown of the DAE
%
%    U = [ u_0 ;  u_C1 ;  u_C2 ; ...  ; u_CH ; u_S1 ;  u_S2 ; ...  ; u_SH ]  Fourier coef (colum) vector for u
%
%    Z= [ P1 ; Q1 ; P2 ; Q2 ; R ; W ; V ]  matrix of unknown Fourier coefficients
%
%    Z= reshape(Z, nz_tot*nb_coef,1)  -> column vector of unknown Fourier coefficients
%
%   If the solution is Periodic, nb_coef = 2*H+1,
%   if the solution is Quasi-periodic, nb_coef = (2*H+2)*H+1

global U Section Diagram    % Global variables to export point from the diagram.

addpath(genpath('../../SRC'))

%% Parameters of the system
nz= 4;     % number of main equations of the DAE system
nz_aux =3; % number of auxiliary equations of the DAE system
H = [5 5]; % number of harmonics
% if H=0, it continues the equilibrium
% if H is an integer, it continues a limit cycle.
% if H is a couple of integer, with H(1)=H(2), it continues a
% quasi-periodic orbit.

% specific parameters
parameters.w1 = 1;
parameters.w2 = 2.5;

parameters.a1=1e-02;
parameters.a2=5e-01;
parameters.a3=2;
parameters.a4=2;

parameters.b1=2.5*parameters.a1;
parameters.b2=2*parameters.a2;
parameters.b3=2*parameters.a3;
parameters.b4=2*parameters.a4;

%% initialization of the system
sys=SystODE(nz,nz_aux,H,@equations,@point_display,@global_display,parameters);

%% starting point

if numel(H) == 2
    %%% Quasi-periodic point :
    %-- Be careful of the displayed variables /!\
    %load('U_H4_QP');
    %%% From Neimarck-Sacker bifurcation
    load('Section_NS1_H10');
    U = sys.init_NS(Section);
    
    U0 = sys.init_Hdiff(U);
    
    % Display the cosine (1,0) and the cosine (0,1) of the third variable with
    % respect to the bifurcation parameter lambda.
    dispvars = [sys.getcoord('lambda') sys.getcoord('sin',3,[1 0]) ; ...
        sys.getcoord('lambda') sys.getcoord('cos',3,[0 1])];
    
    stabcheck = 0;
    
elseif H > 0
        %%% Periodic point :
        omega = 2.5; 
        lambda=0.026; % Behind Hopf bifurcations !
        Z0  =randn(2*H+1,sys.nz_tot)*1e-5;
        Z0(2,1) = .02;
        Z0(:,2) = omega*sys.D(Z0(:,1));
        Z0(2,3) = .2;
        Z0(:,4) = omega*sys.D(Z0(:,1));
    
        Z0(:,5) = sys.Prod(Z0(:,1)+Z0(:,3),Z0(:,1)+Z0(:,3));
        Z0(:,6) = sys.Prod(Z0(:,1)+Z0(:,3),Z0(:,2)+Z0(:,4));
        Z0(:,7) = sys.Prod(Z0(:,5),Z0(:,2)+Z0(:,4));
    
        U0 = sys.init_U0(Z0,omega,lambda);
    
    %%% Initialization after the first Hopf bifurcation
    load('Section_Hopf1');
    U0 = sys.init_Hopf(Section);
%     
    % Display the first cosine and the first sine of the first variable with
    % respect to the bifurcation parameter lambda.
    dispvars = [sys.getcoord('lambda') sys.getcoord('cos',1,1) ; ...
        sys.getcoord('lambda') sys.getcoord('sin',1,1)];
    
    stabcheck = 1;
    
elseif H == 0
    %%% Equilibrium :
    p1 = 0;
    q1 = 0;
    p2 = 0;
    q2 = 0;
    lambda = 0;
    r = 0;
    w = 0;
    v = 0;
    
    U0 = [p1;q1;p2;q2;lambda;r;w;v];
    
    % Display the first and the third variable with respect to the
    % bifurcation parameter lambda.
    dispvars = [sys.neq+1 1 ; sys.neq+1 2];
    
    stabcheck = 1;
    
end


%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'ANMthreshold'    ,1e-8 , ...
    'NRthreshold'     ,1e-5, ...
    'NRstart'         ,0, ...
    'StabilityCheck'  ,stabcheck, ...
    'PointDisplay'    ,1, ...
    'GlobalDisplay'   ,0, ...
    'displayvariables',dispvars);     % MANLAB run


warndlg('Advanced example. Perform 5 continuation steps on the branch before allowing NRcorrections of type 2 (variable parameter). Then the quasi-periodic solution branch will be continued accurately.');
