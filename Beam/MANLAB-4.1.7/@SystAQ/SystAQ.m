classdef SystAQ < Syst
    
    properties
        
        iC=[] ; vC=[] ;                   % list for the sparse (order 1) tensor C
        iL=[] ; jL=[] ; vL=[] ;           % list for the sparse (order 2) tensor L
        iQ=[] ; jQ=[] ; kQ=[] ; vQ=[] ;   % list for the sparse (order 3) tensor Q
        idL=[]; jdL=[]; vdL=[];           % list for the sparse (order 2) tensor dL
        idQ=[]; jdQ=[]; kdQ=[]; vdQ=[];   % list for the sparse (order 3) tensor dQ
        
        writing;
    end
    
    methods
        
        function sys = SystAQ(nz,nz_aux,equations,point_display,global_display,parameters,writing,operators)
            
            if nargin<7; writing = 'standard'; end
            
            neq = nz;
            neq_aux = nz_aux;
            sys=sys@Syst('neq',neq,'neq_aux',neq_aux);
            
            sys.type = 'AQ';
            sys.writing = writing;
            
            sys.parameters = parameters;
            sys.equations = equations;
            sys.point_display = point_display;
            sys.global_display = global_display;
            
            if nargin<8
                switch writing
                    case 'standard'
                        sys = get_operators(sys);
                    case 'vectorial'
                        disp('vectorial initialization of the system.');
                        sys = get_operators_vec(sys);
                end
            else
                sys.iC = operators.iC;
                sys.vC = operators.vC;
                sys.iL = operators.iL;
                sys.jL = operators.jL;
                sys.vL = operators.vL;
                sys.iQ = operators.iQ;
                sys.jQ = operators.jQ;
                sys.kQ = operators.kQ;
                sys.vQ = operators.vQ;
                sys.idL = operators.idL;
                sys.jdL = operators.jdL;
                sys.vdL = operators.vdL;
                sys.idQ = operators.idQ;
                sys.jdQ = operators.jdQ;
                sys.kdQ = operators.kdQ;
                sys.vdQ = operators.vdQ;
            end
            
            % Arclength only on the main variables
            sys.arclengthdef = sparse((1:neq+1),ones(1,neq+1),1,sys.ninc,1);
            
            sys.R = @R;
            
            function [Rf] = R(sys,Uf)
                
                if isempty(sys.idL)
                    Rf = sparse(sys.iC,ones(size(sys.iC)),sys.vC,sys.neq_tot,1) + ...
                        sparse(sys.iL,ones(size(sys.iL)),sys.vL.*Uf(sys.jL),sys.neq_tot,1) + ...
                        sparse(sys.iQ,ones(size(sys.iQ)),sys.vQ.*(Uf(sys.jQ).*Uf(sys.kQ)),sys.neq_tot,1);
                elseif nargin(sys.equations) == 2
                    Rf = sys.equations(sys,Uf);
                elseif nargin(sys.equations) == 3
                    Rf = sys.equations(sys,Uf,zeros(size(Uf)));
                else
                    warndlg('Unsupported file equations.m');
                end
                
            end
            
        end
        
    end
end



