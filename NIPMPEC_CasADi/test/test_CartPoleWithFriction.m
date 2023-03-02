clear all
clc
delete Gen_InitialGuess.mat

%% problem formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 400;
x_Dim = 4;
tau_Dim = 1;
p_Dim = 1;
InitState = [1; 0/180*pi; 0; 0];
RefState = [1; 180/180*pi; 0; 0];
EndState = [1; 180/180*pi; 0; 0];

x = SX.sym('x', x_Dim, 1);
tau = SX.sym('tau', tau_Dim, 1);
p = SX.sym('p', p_Dim, 1); 

% dynamics
mass = [1; 0.1];
linkLength = 1;
g = 9.8;

M = [mass(1) + mass(2),                 mass(2) * linkLength * cos(x(2));...
     mass(2) * linkLength * cos(x(2)),  mass(2) * linkLength^2];
C = [0,   -mass(2) * linkLength * x(4) * sin(x(2));...
     0,   0]; 
G = [0;...
     -mass(2) * g * linkLength * sin(x(2))];
Bu = [tau(1);...
      0]; 
P = [p(1);...
     0]; % friction bewteen cart and ground  
H = G + Bu + P - C * [x(3); x(4)];
f = [x(3:4);...
    inv(M)*H];% xDot = f(x, tau, p)
f_Fun = Function('f_Fun', {x, tau, p}, {f}, {'x', 'tau', 'p'}, {'f'});

% equilibrium constraint
eqlbm.l = -2;
eqlbm.u = 2;
eqlbm.K = x(3);
K_Fun = Function('K_Fun', {x, tau, p}, {eqlbm.K}, {'x', 'tau', 'p'}, {'K'});

% cost function
xWeight_stage = [1; 100; 1; 1];
tauWeight = 1;
L_stageCost = 0.5*(x - RefState)'*diag(xWeight_stage)*(x - RefState) + 0.5*tau'*diag(tauWeight)*tau;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, tau, p}, {L_stageCost}, {'x', 'tau', 'p'}, {'L_stageCost'});

xWeight_terminal = [1; 100; 10; 20];
L_terminalCost = 0.5*(x - EndState)'*diag(xWeight_terminal)*(x - EndState);
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, tau, p}, {L_terminalCost}, {'x', 'tau', 'p'}, {'L_terminalCost'});

% inequality constraint
x_Max = [5; 240/180*pi; 20; 20];
x_Min = [0; -240/180*pi; -20; -20];
tau_Max = 30;
tau_Min = -30;
G_formula =...
    [x_Max - x;...
    x - x_Min;...
    tau_Max - tau;...
    tau - tau_Min];
G_formula_Fun = Function('G_formula_Fun', {x, tau, p}, {G_formula}, {'x', 'tau', 'p'}, {'G'});

% equality constraint
C_formula = [];

% formulate MPEC
X = SX.sym('X', x_Dim, nStages); % state variable
TAU = SX.sym('TAU', tau_Dim, nStages); % control variable
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
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end    
    x_n = X(:, n);
    tau_n = TAU(:, n);
    p_n = P(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, tau_n, p_n);
    L = L + L_n*timeStep;
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, tau_n, p_n);
        L = L + L_terminal;
    end    
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, tau_n, p_n);
    % equality constraint
    C(:, n) = F_n;    
    % inequality constraint
    G_n = G_formula_Fun(x_n, tau_n, p_n);
    G(:, n) = G_n;    
    % equilibrium constraint
    K_n = K_Fun(x_n, tau_n, p_n);
    K(:, n) = K_n;    
end

MPEC.x = reshape([X; TAU], (x_Dim + tau_Dim) * nStages, 1);
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
%%
solver.generateInitialGuess();
Gen_InitialGuess = load('Gen_InitialGuess.mat');
Var_Init = Gen_InitialGuess.Var;

%% solving MPEC
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-6;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-8;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-8;

solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;
solver.Option.linearSystemSolver = 'mldivide_sparse'; % 'linsolve_Sym_dense', 'mldivide_dense', 'mldivide_sparse', 'pinv'

solver.Option.employSecondOrderCorrection = false;
solver.Option.LineSearch.stepSize_Min = 0.01;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-8;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-8;

tic
[solution, Info] = solver.solveMPEC(Var_Init);
toc
%% show result
% iteration process information
solver.showResult(Info)

XTAU_Opt = reshape(solution.x, (x_Dim + tau_Dim), nStages);
x_Opt = XTAU_Opt(1 : x_Dim, :);
tau_Opt = XTAU_Opt(x_Dim + 1 : end, :);
p_Opt = reshape(solution.p, p_Dim, nStages);
timeAxis = 0 : timeStep : nStages * timeStep;

K_Fun_map = K_Fun.map(nStages);
K_value = K_Fun_map(x_Opt, tau_Opt, p_Opt);
K_value = full(K_value);

figure(111)
subplot(4,1,1)
plot(timeAxis, [InitState(1), x_Opt(1, :)], 'r',...
     timeAxis, [InitState(2), x_Opt(2, :)], 'g','LineWidth',1.2)
legend('cart', 'pole')
xlabel('time [s]')
title('position')

subplot(4,1,2)
plot(timeAxis, [InitState(3), x_Opt(3, :)], 'r',...
     timeAxis, [InitState(4), x_Opt(4, :)], 'g', 'LineWidth',1.2)
xlabel('time [s]')
title('velocity')

subplot(4,1,3)
plot(timeAxis(2:end), tau_Opt(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('control')

subplot(4,1,4)
plot(timeAxis(2:end), p_Opt(1, :), 'k',...
     timeAxis(2:end), K_value(1, :), 'b', 'LineWidth', 1.2)
legend('friction', 'cart vel') 
xlabel('time [s]')
title('equilibrium dynamics')