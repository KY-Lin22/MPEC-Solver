% example from paper: D.Fernandez et al, Stabilized sequential quadratic
% programming for optimization and a stabilized Newton-type method for
% variational problems, Mathematical Programming, 2010
clear all
clc
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

% degenerate NLP problem formulation
NLP = struct;
x = SX.sym('x', 3, 1); % optimal variable x
NLP.x = x;
NLP.f = x(1)*x(2) - 1/2*x(2)^2; % cost function f(x)
NLP.h = x(3); % equality constraint h(x) = 0 (h can not be empty in current solver version because 1 x 1 + 0 x 0 = 0 x 0 in Lagrangian)
NLP.g = [x(2)^2;...
    -2*x(1) + x(2);...
    x(1)- 2*x(2)]; % inequality constraint g(x) <= 0

% provide initial guess for x
x_Init = [10; 10];

% create stabilized_SQP_Izamailov2015 object
solver = stabilized_SQP_Izamailov2015(NLP);

% solve degenerate NLP
[x_Opt, Info] = solver.solveNLP([x_Init; 0]);

disp('x_Opt: ')
disp(x_Opt(1:2))
disp('mu: ')
disp(Info.mu)
