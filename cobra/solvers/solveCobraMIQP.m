function solution = solveCobraMIQP(MIQPproblem,varargin)
%solveCobraQP Solve constraint-based QP problems
%
% solution = solveCobraQP(MIQPproblem,solver,verbFlag,solverParams)
%
% % Solves problems of the type 
%
%      min   osense * 0.5 x' * F * x + osense * c' * x
%      s/t   lb <= x <= ub
%            A * x  <=/=/>= b
%            xi = integer
%
%INPUT
%MIQPproblem    Structure containing the following fields describing the QP
%               problem to be solved
%  A                LHS matrix
%  b                RHS vector
%  F                F matrix for quadratic objective (see above)
%  c                Objective coeff vector
%  lb               Lower bound vector
%  ub               Upper bound vector
%  osense           Objective sense (-1 max, +1 min)
%  csense           Constraint senses, a string containting the constraint 
%                   sense for each row in A ('E', equality, 'G' greater 
%                   than, 'L' less than).
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
% the solver defined in the CBT_MIQP_SOLVER global variable (set using
% changeCobraSolver). Solvers currently available are 'tomlab_cplex'
%
%OUTPUT
% solution      Structure containing the following fields describing a QP
%               solution
%  full             Full QP solution vector
%  obj              Objective value
%  solver           Solver used to solve QP problem
%  stat             Solver status in standardized form (see below)
%                       1   Optimal solution found
%                       2   Unbounded solution
%                       0   Infeasible QP
%                      -1   No optimal solution found (time limit etc)
%                       3   Solution exists but with problems
%  origStat         Original status returned by the specific solver 
%  time             Solve time in seconds
%
%
% Markus Herrgard 6/8/07

global CBT_MIQP_SOLVER;
solver = CBT_MIQP_SOLVER;

%optional parameters
optParamNames = {'printLevel', 'saveInput', 'timeLimit'};
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
[printLevel, saveInput, timeLimit] = getCobraSolverParams('QP',optParamNames,parameters);

[A,b,F,c,lb,ub,csense,osense, vartype] = ...
    deal(MIQPproblem.A,MIQPproblem.b,MIQPproblem.F,MIQPproblem.c,MIQPproblem.lb,MIQPproblem.ub,...
    MIQPproblem.csense,MIQPproblem.osense, MIQPproblem.vartype);

t_start = clock;
switch solver
%% CPLEX through TOMLAB
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
        intVars = find((vartype == 'B') | (vartype == 'I'));
        %tomlabProblem = qpAssign(osense*F,osense*c,A,b_L,b_U,lb,ub,[],'CobraQP');
        tomlabProblem  = miqpAssign(osense*F, osense*c, A, b_L, b_U, lb, ub,[], ...
                             intVars, [],[],[],'CobraMIQP');
        tomlabProblem.CPLEX.LogFile = 'MIQPproblem.log';
        
        %optional parameters
        PriLvl = printLevel;
        
        %Save Input if selected
        if ~isempty(saveInput)
            fileName = parameters.saveInput;
            if ~find(regexp(fileName,'.mat'))
                fileName = [fileName '.mat'];
            end
            display(['Saving MIQPproblem in ' fileName]);
            save(fileName,'MIQPproblem')
        end
        tomlabProblem.MIP.cpxControl.TILIM = timeLimit; % time limit
        tomlabProblem.MIP.cpxControl.THREADS = 1; % by default use only one thread
        Result = tomRun('cplex', tomlabProblem, PriLvl);
        
        x = Result.x_k;
        f = osense*Result.f_k;
        stat = Result.Inform;
        if (stat == 1 ||stat == 101 || stat == 102)
            solStat = 1; % Optimal
        elseif (stat == 3 || stat == 4)
            solStat = 0; % Infeasible
        elseif (stat == 103)
            solStat = 0; % Integer Infeasible
        elseif (stat == 2 || stat == 118 || stat == 119)
            solStat = 2; % Unbounded
        elseif (stat == 106 || stat == 108 || stat == 110 || stat == 112 || stat == 114 || stat == 117)
            solStat = -1; % No integer solution exists
        elseif (stat >= 10)
            solStat = -1; % No optimal solution found (time or other limits reached, other infeasibility problems)
        else
            solStat = 3; % Solution exists, but either scaling problems or not proven to be optimal
        end 
    otherwise
        error(['Unknown solver: ' solver]);
end
%%
t = etime(clock, t_start);

solution.obj = f;
solution.solver = solver;
solution.stat = solStat;
solution.origStat = stat; 
solution.time = t;
solution.full = x;