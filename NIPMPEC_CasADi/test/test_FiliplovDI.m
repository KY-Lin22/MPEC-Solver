clear all
clc
delete Gen_InitialGuess.mat

%% Problem Formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
% FilippovDI: an example from Stewart's paper:
% Optimal control of systems with discontinuous differential equations


timeStep = 0.02;
nStages = 100;
x_Dim = 1;
p_Dim = 1;
InitState = -1;
RefState = 0;
EndState = 5/3;

x = SX.sym('x', x_Dim, 1);
p = SX.sym('p', p_Dim, 1); % variable y in the Stewart's paper

% dynamics
f1 = 1;% switch function > 0
f2 = 3; % switch function < 0
f = f1*(1 - p) + f2*p;
f_Fun = Function('f_Fun', {x, p}, {f}, {'x', 'p'}, {'f'});

% equilibrium constraint
eqlbm.l = 0;
eqlbm.u = 1;
eqlbm.K = x;
K_Fun = Function('K_Fun', {x, p}, {eqlbm.K}, {'x', 'p'}, {'K'});

% cost function
L_stageCost = (x - RefState)^2;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, p}, {L_stageCost}, {'x', 'p'}, {'L_stageCost'});

L_terminalCost = (x - EndState)^2;
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, p}, {L_terminalCost}, {'x', 'p'}, {'L_terminalCost'});

% inequality constraint
x_Max = 10000;
x_Min = -10000;
G_formula =...
    [x_Max - x;...
    x - x_Min];
G_formula_Fun = Function('G_formula_Fun', {x, p}, {G_formula}, {'x', 'p'}, {'G'});

% equality constraint
C_formula = [];

% formulate MPEC
X = SX.sym('X', x_Dim, nStages); % optimal variable
P = SX.sym('P', p_Dim, nStages); % algebraic variable
L = 0; % init cost function
G_Dim = size(G_formula, 1);
G = SX.sym('G', G_Dim, nStages); % init inequality constraint
C_Dim = size([C_formula; f], 1);
C = SX.sym('C', C_Dim, nStages); % init equality constraint
K_Dim = size(eqlbm.K, 1);
K = SX.sym('K', K_Dim, nStages); % init function K
l = repmat(eqlbm.l, 1, nStages);
u = repmat(eqlbm.u, 1, nStages);
for n = 1 : nStages
    % variable
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end
    x_n = X(:, n);
    p_n = P(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, p_n);
    L = L + L_n*timeStep;
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, p_n);
        L = L + L_terminal;
    end
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, p_n);
    % equality constraint
    C(:, n) = F_n;
    % inequality constraint
    G_n = G_formula_Fun(x_n, p_n);
    G(:, n) = G_n;
    % equilibrium constraint
    K_n = K_Fun(x_n, p_n);
    K(:, n) = K_n;
end

MPEC.x = reshape(X, x_Dim * nStages, 1);
MPEC.p = reshape(P, p_Dim * nStages, 1);
MPEC.L = L;
MPEC.G = reshape(G, G_Dim * nStages, 1);
MPEC.C = reshape(C, C_Dim * nStages, 1);
MPEC.l = reshape(l, [], 1);
MPEC.u = reshape(u, [], 1);
MPEC.K = reshape(K, K_Dim * nStages, 1);

%% create Solver
% create solver object
solver = NIPMPEC_CasADi(MPEC);

solver.showInfo();

solver.generateInitialGuess();
Gen_InitialGuess = load('Gen_InitialGuess.mat');
Var_Init = Gen_InitialGuess.Var;

%% solving MPEC
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-4;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-6;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-6;

solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;
solver.Option.linearSystemSolver = 'mldivide_sparse'; % 'linsolve_Sym_dense', 'mldivide_dense', 'mldivide_sparse', 'pinv'

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-3;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-3;

tic
[solution, Info] = solver.solveMPEC(Var_Init);
toc

%% show result
% iteration process information
solver.showResult(Info)
% plot solution
timeAxis = 0 : timeStep : nStages * timeStep;

f_Fun_map = f_Fun.map(nStages);
f_value = f_Fun_map(solution.x', solution.p');
f_value = full(f_value);

figure(112)
subplot(3,2,1)
plot(timeAxis, [InitState, solution.x'], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time(s)')
title('system state')

subplot(3,2,2)
plot(timeAxis(2:end), f_value, 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(3,2,3)
plot(timeAxis(2:end), solution.p', 'k', 'LineWidth', 1.2)
legend('y') 
xlabel('time(s)')
title('smoothing function')

subplot(3,2,4)
plot(timeAxis(2:end), solution.x', 'g', 'LineWidth', 1.2)
legend('z := x') 
xlabel('time(s)')
title('switch function')

subplot(3,2,5)
plot(timeAxis(2:end), solution.x', 'k',...
    timeAxis(2:end), solution.p', 'b', 'LineWidth', 1.2)
legend('z(K)', 'y(p)') 
xlabel('time(s)')
title('checking BVI')
