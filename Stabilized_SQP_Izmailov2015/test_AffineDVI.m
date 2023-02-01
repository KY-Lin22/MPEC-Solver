% Affine DVI

clear all
clc
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.1;
nStages = 3; 
s = 1e-1; % slack 

%% Dynamics
% dynamics variables
tau_Dim = 1;% control dim
x_Dim = 2; % state dim
p_Dim = 1; % equilibrium var dim
w_Dim = p_Dim; % auxiliary var dim
u_Dim = tau_Dim + p_Dim + w_Dim; % u = [tau; p; w]

x = SX.sym('x', x_Dim);
u = SX.sym('u', u_Dim);
tau = u(1 : tau_Dim);
p = u(tau_Dim + 1 : tau_Dim + p_Dim);
w = u(tau_Dim + p_Dim + 1 : end);

% state equations
A = [1, -3; ...
    -8, 10];
B = [-3;...
    -1];
F = [4;...
    8];
f = A * x + B * p + F * tau; % xDot = f(tau, x, p)
f_func = Function('f_func',{x,u}, {f}, {'x', 'u'}, {'f'});

% equilibrium dynamics
lbp = -1;
ubp = 1;
C = [1, -3];
D = 5;
E = 3;
K = C * x + D * p + E * tau;
K_func = Function('K_func',{x,u}, {K}, {'x', 'u'}, {'K'});

% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
BVI_eq = w - K;% auxiliary variable, BVI_eq = 0
BVI_ineq = [lbp - p;...
    p - ubp;...%  p in [l, u] 
    (p - lbp) * w - s;...
    -(ubp - p) * w - s]; % regularization (Note: BVI_ineq <= 0)

BVI_func = Function('BVI_func', {x,u}, {BVI_eq, BVI_ineq}, {'x', 'u'}, {'BVI_eq', 'BVI_ineq'});

% optimal variable bounds
tau_Max = 5;
tau_Min = -5;
x_Max = [5; 5];
x_Min = [-5; -5];

optVarBounds = [x - x_Max;...
    x_Min - x;...
    tau - tau_Max;...
    tau_Min - tau]; % optVarBounds <= 0
optVarBounds_func = Function('optVarBounds_func', {x, u}, {optVarBounds}, {'x', 'u'}, {'optVarBounds'});

%% OCPEC
% specify initial and end state, cost ref and weight matrix
InitState = [-1/2; -1];
EndState = [0; 0];

StageCost.xRef = repmat(EndState, 1, nStages);
StageCost.tauRef = zeros(1, nStages);
StageCost.xWeight = [20; 20];
StageCost.tauWeight = 1;
TerminalCost.xRef = EndState;
TerminalCost.tauRef = 0;
TerminalCost.xWeight = [20; 20];
TerminalCost.tauWeight = 1;

pWeight = 0.001 * eye(p_Dim);
wWeight = 0.001 * eye(w_Dim);

% optimal variables
X = SX.sym('X', x_Dim, nStages);
U = SX.sym('U', u_Dim, nStages); 

% cost function and constraints
L = 0; % initial cost function
h_Dim = size([f; BVI_eq], 1);
g_Dim = size([BVI_ineq; optVarBounds], 1);
h = SX.sym('h', h_Dim, nStages);
g = SX.sym('g', g_Dim, nStages);

for n = 1 : nStages
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end
    x_n = X(:, n);
    u_n = U(:, n);
    tau_n = u_n(1 : tau_Dim, 1);
    p_n = u_n(tau_Dim + 1 : tau_Dim + p_Dim, 1);
    w_n = u_n(tau_Dim + p_Dim + 1 : end, 1);    
    % cost function
    L_n = 0.5 * (x_n - StageCost.xRef(:, n))' * diag(StageCost.xWeight) * (x_n - StageCost.xRef(:, n))...
        + 0.5 * (tau_n - StageCost.tauRef(:, n))' * diag(StageCost.tauWeight) * (tau_n - StageCost.tauRef(:, n))...
        + 0.5 * p_n' * pWeight * p_n...
        + 0.5 * w_n' * wWeight * w_n;
    L = L + L_n * timeStep;
    if n == nStages
        L_Terminal = 0.5 * (x_n - TerminalCost.xRef)' * diag(TerminalCost.xWeight) * (x_n - TerminalCost.xRef)...
            + 0.5 * (tau_n - TerminalCost.tauRef)' * diag(TerminalCost.tauWeight) * (tau_n - TerminalCost.tauRef);
        L = L + L_Terminal;
    end    
    % discretize dynamics by implicit euler method
    F_n = x_nPrev - x_n + timeStep * f_func(x_n, u_n);
    % reformulated equilibrium constraint
    [BVI_eq_n, BVI_ineq_n] = BVI_func(x_n, u_n);
    % optimal variable bounds
    optVarBounds_n = optVarBounds_func(x_n, u_n);
    % constraint function
    h(:, n) = [F_n;...
        BVI_eq_n];
    g(:, n) = [BVI_ineq_n;...
        optVarBounds_n];
end

%% Solver
% reshape optimal variable and constraint
XU = reshape([X;U], (x_Dim + u_Dim) * nStages, 1);
h = reshape(h, h_Dim * nStages, 1);
g = reshape(g, g_Dim * nStages, 1);

NLP = struct;
NLP.x = XU;
NLP.f = L; % cost function f(x)
NLP.h = h; % equality constraint h(x) = 0 (h can not be empty in current solver version because 1 x 1 + 0 x 0 = 0 x 0 in Lagrangian)
NLP.g = g; % inequality constraint g(x) <= 0

% provide initial guess for x
x_Init = zeros((x_Dim + u_Dim) * nStages, 1);

% create stabilized_SQP_Izamailov2015 object
solver = stabilized_SQP_Izamailov2015(NLP);

% solve degenerate NLP
[x_Opt, Info] = solver.solveNLP(x_Init);

x_Opt = reshape(full(x_Opt), x_Dim + u_Dim, nStages);
figure(1)
plot([InitState(1), x_Opt(1, :)], '-g')
hold on
plot([InitState(2), x_Opt(2, :)], '-r')
hold on
plot(x_Opt(3, :), '-k')
legend('x1', 'x2', 'tau')
