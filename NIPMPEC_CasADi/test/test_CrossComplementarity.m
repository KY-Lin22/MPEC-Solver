clear all
clc
delete Gen_InitialGuess.mat
%%
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
%
% min (x_1 - 1)^2 + (x_2 - 1)^2 + (x_3 - 1)^2 + (x_4 - 1)^2
% s.t. 0 <= x_1 >= x_2 >= 0
%      0 <= x_3 >= x_4 >= 0
%      x_1x_4 = 0
%      x_2x_3 = 0

% symbolic respresentation of variable
x = SX.sym('x', 2, 1); % x(1)denotes x_1; 
                       % x(2)denotes x_2
p = SX.sym('p', 2, 1); % p(1)denotes x_3; 
                       % p(2)denotes x_4
                       
% cost function L, inequality constraint G >= 0 and equality constrain C = 0 
L = (x(1) - 1)^2 + (x(2) - 1)^2 + (p(1) - 1)^2 + (p(2) - 1)^2;
G = [x(1) - x(2);...% denotes x_1 >= x_2 in cross complementary constraint
    p(1) - p(2)]; % denotes x_3 >= x_4 in cross complementary constraint
C = [];

% equilibrium constraint (element-wise)
% NOTE: if the equilibrium constraint is nonlinear complementary constraint(NCP, p>=0, K>=0, pK = 0, p is the algebraic variable, K is the function),
%       then you can formulate NCP by step 1: providing the symbolic respresentation of p and K,
%                                     step 2: setting l = 0 and u = inf,
%       and you DO NOT need to formulate the inequality constraint p >= 0 and K >= 0 (which will be automatically formulated by solver)
l = [0;...
    0];% specify l = 0, u = inf to formulate two sets of NCP: p(1) and K(1), p(2) and K(2)
u = [Inf;...
    Inf];
K = [x(2);...% K(1) is complementary to p(1), so let K(1):= x(2) denoting x_2x_3 = 0
    x(1)];   % K(2) is complementary to p(2), so let K(2):= x(1) denoting x_1x_4 = 0

% create problem struct
MPEC.x = x;
MPEC.p = p;
MPEC.L = L;
MPEC.G = G;
MPEC.C = C;
MPEC.l = l;
MPEC.u = u;
MPEC.K = K;

% create solver object
solver = NIPMPEC_CasADi(MPEC);

solver.showInfo();
%% generate and load initial guess
solver.generateInitialGuess();
Gen_InitialGuess = load('Gen_InitialGuess.mat');
Var_Init = Gen_InitialGuess.Var;

%% solve
% solver option
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-6;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-8;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-8;

solver.Option.employSparsePattern = true;
solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-8;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-8;

tic
[solution, Info] = solver.solveMPEC(Var_Init);
toc

x_1 = solution.x(1)
x_2 = solution.x(2)
x_3 = solution.p(1)
x_4 = solution.p(2)
% solver.showResult(Info)