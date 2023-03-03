clear all
clc
delete Gen_InitialGuess.mat
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

x = SX.sym('x', 1, 1);
p = SX.sym('p', 1, 1);
L = (x(1) - 1)^2 + (p(1) - 1)^2;
G = [];
C = [];
l = 0;
u = Inf;
K = x(1);

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

solver.generateInitialGuess();
Gen_InitialGuess = load('Gen_InitialGuess.mat');
Var_Init = Gen_InitialGuess.Var;

%%
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-4;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-6;

solver.Option.employSparsePattern = true;
solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;

solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-3;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-3;

tic
[solution, Info] = solver.solveMPEC(Var_Init);
toc
% solver.showResult(Info)