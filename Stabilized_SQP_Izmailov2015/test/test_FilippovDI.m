% Filippov DI

clear all
clc
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.02;
nStages = 100; 
s = 0; % slack 

%% Dynamics
x_Dim = 1;
u_Dim = 1; % y

x = SX.sym('x', x_Dim);
u = SX.sym('u', u_Dim);

% state equations
f1 = 1; % switch function > 0
f2 = 3; % switch function < 0
f = f1 * (1 - u) + f2 * u;
f_func = Function('f_func',{x,u}, {f}, {'x', 'u'}, {'f'});

% equilibrium dynamics (Note: BVI_ineq <= 0)
lbp = 0;
ubp = 1;
BVI_ineq = [lbp - u;...
    u - ubp;...%  y in [l, u] 
    (u - lbp) * x - s;...
    -(ubp - u) * x - s]; % regularization 
BVI_func = Function('BVI_func', {x,u}, {BVI_ineq}, {'x', 'u'}, {'BVI_ineq'});

%% OCPEC
% specify initial, cost ref and weight matrix
InitState = -1;

StageCost.xWeight = 2;
TerminalCost.xRef = 5/3;
TerminalCost.xWeight = 2;

% optimal variables
X = SX.sym('X', x_Dim, nStages);
U = SX.sym('U', u_Dim, nStages); 

% cost function and constraints
L = 0; % initial cost function
h_Dim = size(f, 1);
g_Dim = size(BVI_ineq, 1);
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
    % cost function
    L_n = 0.5 * (x_n)' * diag(StageCost.xWeight) * (x_n);
    L = L + L_n * timeStep;
    if n == nStages
        L_Terminal = 0.5 * (x_n - TerminalCost.xRef)' * diag(TerminalCost.xWeight) * (x_n - TerminalCost.xRef);
        L = L + L_Terminal;
    end    
    % discretize dynamics by implicit euler method
    F_n = x_nPrev - x_n + timeStep * f_func(x_n, u_n);
    % reformulated equilibrium constraint
    BVI_ineq_n = BVI_func(x_n, u_n);
    % constraint function
    h(:, n) = F_n;
    g(:, n) = BVI_ineq_n;
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
plot([x_Opt(2, :)], '-r')
legend('x', 'u')