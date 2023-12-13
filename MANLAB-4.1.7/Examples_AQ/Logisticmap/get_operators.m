function [operators] = get_operators(parameters)

operators.iC=[];
operators.vC=[];
operators.idL=[];
operators.jdL=[];
operators.vdL=[];
operators.idQ=[];
operators.jdQ=[];
operators.kdQ=[];
operators.vdQ=[];

N = parameters.N;

xind = (1:N-1)';
xp1 = (2:N)';

operators.iC = [1 ; N+1];
operators.vC = [-parameters.x0;-1];

operators.iL = [1; xp1];
operators.jL = [1 ; xind];
operators.vL = [1;-ones(N-1,1)];

operators.iQ = [xp1;              xp1;xp1 ; N+1 ; N+1];
operators.jQ = [xp1;              xind;ones(N-1,1)*(N+2); N+1 ; N+2];
operators.kQ = [ones(N-1,1)*(N+2);xind;xp1; N+2 ; N+1];
operators.vQ = [0.5*ones(N-1,1);ones(N-1,1);0.5*ones(N-1,1) ;.5 ; .5];

end

