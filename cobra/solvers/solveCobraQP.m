function solution = solveCobraQP(QPproblem,varargin)
%solveCobraQP Solve constraint-based QP problems
%
% solution = solveCobraQP(QPproblem,parameters)
%
% % Solves problems of the type 
%
%      min   0.5 x' * F * x + osense * c' * x
%      s/t   lb <= x <= ub
%            A * x  <=/=/>= b
%
%INPUT
% QPproblem Structure containing the following fields describing the QP
% problem to be solved
%  A      LHS matrix
%  b      RHS vector
%  F      F matrix for quadratic objective (must be positive definite)
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).
%
%OPTIONAL INPUTS
% Optional parameters can be entered using parameters structure or as
% parameter followed by parameter value: i.e. ,'printLevel',3)
%
% parameters    Structure containing optional parameters as fields.
%  printLevel   Print level for solver
%  saveInput    Saves LPproblem to filename specified in field. 
%               Setting parameters = 'default' uses default setting set in
%               getCobraSolverParameters.
%
% The solver is defined in the CBT_MILP_SOLVER global variable
% (set using changeCobraSolver). Solvers currently available are
% 'tomlab_cplex', 'mosek' and 'qpng' (limited support for small problems)
%
%OUTPUT
% solution  Structure containing the following fields describing a QP
%           solution
%  full     Full QP solution vector
%  obj      Objective value
%  solver   Solver used to solve QP problem
%  stat     Solver status in standardized form (see below)
%  origStat Original status returned by the specific solver 
%  time     Solve time in seconds
%
%  stat     Solver status in standardized form
%           1   Optimal solution
%           2   Unbounded solution
%           0   Infeasible
%           -1  No solution reported (timelimit, numerical problem etc)
%
% Markus Herrgard        6/8/07
% Ronan Fleming         12/07/09  Added support for mosek
% Ronan Fleming         18 Jan 10 Added support for pdco
% Josh Lerman           04/17/10 changed def. parameters, THREADS, QPMETHOD

global CBT_QP_SOLVER;

if (~isempty(CBT_QP_SOLVER))
    solver = CBT_QP_SOLVER;
else
    error('No solver found');
end

optParamNames = {'printLevel','saveInput'};
parameters = '';
if nargin ~=1
    if mod(length(varargin),2)==0
        for i=1:2:length(varargin)-1
            if ismember(varargin{i},optParamNames)
                parameters.(varargin{i}) = varargin{i+1};
            else
                error([varargin{i} ' is not a valid optional parameter']);
            end
        end
    elseif strcmp(varargin{1},'default')
        parameters = 'default';
    elseif isstruct(varargin{1})
        parameters = varargin{1};
    else
        display('Warning: Invalid number of parameters/values')
        solution=[];
        return;
    end
end

% Defaults in case the solver does not return anything
x = [];
y = [];
w = [];
f = [];
xInt = [];
xCont = [];
stat = -99;
solStat = -99;

%parameters
[printLevel, saveInput] = getCobraSolverParams('QP',optParamNames,parameters);

[A,b,F,c,lb,ub,csense,osense] = ...
    deal(QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,...
    QPproblem.csense,QPproblem.osense);

t_start = clock;
switch solver
%%
    case 'tomlab_cplex'
        if (~isempty(csense))
            b_L(csense == 'E') = b(csense == 'E');
            b_U(csense == 'E') = b(csense == 'E');
            b_L(csense == 'G') = b(csense == 'G');
            b_U(csense == 'G') = inf;
            b_L(csense == 'L') = -inf;
            b_U(csense == 'L') = b(csense == 'L');
        else
            b_L = b;
            b_U = b;
        end
        tomlabProblem = qpAssign(F,osense*c,A,b_L,b_U,lb,ub,[],'CobraQP');
        
        %optional parameters
        tomlabProblem.PriLvl=printLevel;
        tomlabProblem.MIP.cpxControl.QPMETHOD = 1;
        tomlabProblem.MIP.cpxControl.THREADS = 1;
                
        %Save Input if selected
        if ~isempty(saveInput)
            fileName = saveInput;
            if ~find(regexp(fileName,'.mat'))
                fileName = [fileName '.mat'];
            end
            display(['Saving QPproblem in ' fileName]);
            save(fileName,'QPproblem')
        end
        
        Result = tomRun('cplex', tomlabProblem);
        x = Result.x_k;
        f = osense*Result.f_k;
        origStat = Result.Inform;
        if (origStat == 1)
            stat = 1; % Optimal
        elseif (origStat == 3 || origStat == 4)
            stat = 0; % Infeasible
        elseif (origStat == 2)
            stat = 2; % Unbounded
        elseif (origStat >= 10)
            stat = -1; % No optimal solution found (time or other limits reached, other infeasibility problems)
        else
            stat = 3; % Solution exists, but either scaling problems or not proven to be optimal
        end
        %%
    case 'cplex_direct'
        %% Tomlab cplex.m direct
        %Used with the current script, only some of the control affoarded with
        %this interface is provided. Primarily, this is to change the print
        %level and whether to minimise the Euclidean Norm of the internal
        %fluxes or not.
        %See solveCobraLPCPLEX.m for more refined control of cplex
        %Ronan Fleming 11/12/2008
               
        solution=solveCobraLPCPLEX(QPproblem,printLevel,[],[],[],minNorm);
        %%
   case 'qpng'
        % qpng.m This file is part of GLPKMEX.
        % Copyright 2006-2007 Nicolo Giorgetti.
        %
        % Solves the general quadratic program
        %      min 0.5 x'*H*x + q'*x
        %       x
        % subject to
        %      A*x [ "=" | "<=" | ">=" ] b
        %      lb <= x <= ub
        ctype=csense;
        ctype(('G'==csense))='L';
        ctype(('E'==csense))='E';
        ctype(('L'==csense))='U';
        
        x0=ones(size(QPproblem.A,2),1);
        %equality constraint matrix must be full row rank
        [x, f, y, info] = qpng (QPproblem.F, QPproblem.c*QPproblem.osense, full(QPproblem.A), QPproblem.b, ctype, QPproblem.lb, QPproblem.ub, x0);
        
        f = 0.5*x'*QPproblem.F*x + c'*x;
        
        w=[];
        
        if (info.status == 0)
            stat = 1;
        elseif (info.status == 1)
            stat = 0;
        else
            stat = -1;
        end
        origStat=info.status;
        %%
    case 'mosek'
        if (~isempty(csense))
            b_L(csense == 'E') = b(csense == 'E');
            b_U(csense == 'E') = b(csense == 'E');
            b_L(csense == 'G') = b(csense == 'G');
            b_U(csense == 'G') = inf;
            b_L(csense == 'L') = -inf;
            b_U(csense == 'L') = b(csense == 'L');
        else
            b_L = b;
            b_U = b;
        end
                
        if printLevel>0
            cmd='minimize';
        else
            cmd='minimize echo(0)';
        end
        
        % Optimize the problem.
        % min 0.5*x'*F*x + osense*c'*x
        % st. blc <= A*x <= buc
        %     bux <= x   <= bux
        [res] = mskqpopt(F,osense*c,A,b_L,b_U,lb,ub,[],cmd);
        
        if isempty(res)
            stat=3;
        else
            if isfield(res,'sol')
                origStat=res.sol.itr.solsta;
                if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE') &&  (strcmp(res.sol.itr.solsta,'OPTIMAL') || strcmp(res.sol.itr.solsta,'NEAR_OPTIMAL'))
                    stat=1;
                    % x solution.
                    x = res.sol.itr.xx;
                    f = 0.5*x'*F*x + c'*x;
                    
                    %dual to equality
                    y=res.sol.itr.y;
                    
                    %dual to lower and upper bounds
                    w=res.sol.itr.slx - res.sol.itr.sux;
                else
                    stat=3;
                end
            else
                stat=3;
                origStat=res.rmsg;
            end
        end
        % stat   Solver status
        %           1   Optimal solution found
        %           2   Unbounded solution
        %           0   Infeasible QP
        %           3   Other problem (time limit etc)
        %%
      case 'pdco'
        %-----------------------------------------------------------------------
        % pdco.m: Primal-Dual Barrier Method for Convex Objectives (16 Dec 2008)
        %-----------------------------------------------------------------------
        % AUTHOR:
        %    Michael Saunders, Systems Optimization Laboratory (SOL),
        %    Stanford University, Stanford, California, USA.
        %Interfaced with Cobra toolbox by Ronan Fleming, 18 Jan 2010
        [nMet,nRxn]=size(A);
        d1=ones(nRxn,1)*1e-4;
        %dont minimise the norm of reactions in linear objective
        d1(c~=0)=0;
        d2=1e-5;
        options = pdcoSet;
        
        x0 = ones(nRxn,1);
        y0 = ones(nMet,1);
        z0 = ones(nRxn,1);
        xsize = 1000;
        zsize = 1000;
        options.Method=2; %QR
        options.MaxIter=100;
        options.Print=printLevel;
        %get handle to helper function for objective
        pdObjHandle = @(x) QPObj(x);
        %solve the QP
        [x,y,w,inform,PDitns,CGitns,time] = ...
            pdco(pdObjHandle,A,b,lb,ub,d1,d2,options,x0,y0,z0,xsize,zsize);
        f= c'*x + 0.5*x'*F*x;
        % inform = 0 if a solution is found;
%        = 1 if too many iterations were required;
%        = 2 if the linesearch failed too often;
%        = 3 if the step lengths became too small;
%        = 4 if Cholesky said ADDA was not positive definite.
        if (inform == 0)
            stat = 1;
        elseif (inform == 1 || inform == 2 || inform == 3)
            stat = 0;
        else
            stat = -1;
        end
        origStat=inform;
    %%
    otherwise
        error(['Unknown solver: ' solver]);
end
%%
t = etime(clock, t_start);
solution.obj = f;
solution.solver = solver;
solution.stat = stat;
solution.origStat = origStat; 
solution.time = t;
solution.full = x;
solution.dual = y;
solution.rcost = w;

%Helper function for pdco
%%
    function [obj,grad,hess] = QPObj(x)
        obj  = c'*x + 0.5*x'*F*x;
        grad = c + F*x;
        hess = F;
    end
end