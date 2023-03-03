clc
clear all
%% problem formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.02;
nStages = 100;
sInit = 1e-1;
sEnd = 1e-5; 
zEnd = 1e-5;

x_Dim = 1;
p_Dim = 1;
w_Dim = p_Dim; % auxiliary var
s_Dim = 1; % regularization parameter
x = SX.sym('x', x_Dim);
p = SX.sym('p', p_Dim);
w = SX.sym('w', w_Dim);
s = SX.sym('s', s_Dim);

InitState = -1;
RefState = 0;
EndState = 5/3;

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

% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
BVI = [p - eqlbm.l;...
    eqlbm.u - p;...%  p in [l, u]   
    w - eqlbm.K;...% auxiliary variable
    s - (p - eqlbm.l) * w;...
    s + (eqlbm.u - p) * w];% regularization
BVI_Fun = Function('BVI_Fun', {x, p, w, s}, {BVI}, {'x', 'p', 'w', 's'}, {'BVI'});
lbg_BVI = [0; 0; 0; 0; 0];
ubg_BVI = [Inf; Inf; 0; Inf; Inf];

% cost function
L_stageCost = (x - RefState)^2;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, p}, {L_stageCost}, {'x', 'p'}, {'L_stageCost'});

L_terminalCost = (x - EndState)^2;
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, p}, {L_terminalCost}, {'x', 'p'}, {'L_terminalCost'});

% formulate NLP
X = SX.sym('X', x_Dim, nStages); % state variable
P = SX.sym('P', p_Dim, nStages); % algebraic variable
W = SX.sym('W', w_Dim, nStages); % auxilary variable
S = SX.sym('S', s_Dim, 1);% regularizaton parameter

lbx = -Inf * ones(x_Dim + p_Dim + w_Dim, nStages);
ubx = Inf * ones(x_Dim + p_Dim + w_Dim, nStages);

L = 0; % init cost function
g_Dim = size([f; BVI], 1);
g = SX.sym('g', g_Dim, nStages); % constraint function
lbg = zeros(g_Dim, nStages);
ubg = zeros(g_Dim, nStages);
lbg(size(f, 1) + 1 : end, :) = repmat(lbg_BVI, 1, nStages);
ubg(size(f, 1) + 1 : end, :) = repmat(ubg_BVI, 1, nStages);

for n = 1 : nStages
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end    
    x_n = X(:, n);
    p_n = P(:, n);
    w_n = W(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, p_n);
    L = L + L_n*timeStep;    
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, p_n);
        L = L + L_terminal;
    end 
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, p_n);   
    % reformulated equilibrium constraint
    BVI_n = BVI_Fun(x_n, p_n, w_n, S); 
    % constraint function
    g(:, n) = [F_n;...
        BVI_n];     
end

OCPEC = struct('x', reshape([X; P; W], (x_Dim + p_Dim + w_Dim) * nStages, 1),...,...
    'f', L,...  
    'g', reshape(g, g_Dim * nStages, 1),...
    'p', S);

args.lbx = reshape(lbx, [], 1);
args.ubx = reshape(ubx, [], 1);
args.lbg = reshape(lbg, [], 1);
args.ubg = reshape(ubg, [], 1);
% initialize x0 and s
args.x0 = zeros((x_Dim + p_Dim + w_Dim) * nStages, 1);
args.p = sInit;

% option (solver.print_options)
Option = struct;
Option.print_time = false;
Option.ipopt.max_iter = 2000;
Option.ipopt.tol = 1e-4;
Option.ipopt.mu_target = 0.5 * (zEnd)^2;
Option.ipopt.print_level = 0;

% create solver
solver = nlpsol('solver', 'ipopt', OCPEC, Option);

%% homotopy
kappa_s_times = 0.2;
kappa_s_exp = 1.5;
homotopy_counter = 1;

totalTime_Start = tic;
while true
    % Homotopy (outer) iteration
    HomotopyTime_Start = tic;
    solution = solver('x0', args.x0, 'p',args.p,...
        'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg);
    HomotopyTime = toc(HomotopyTime_Start);
    % check IPOPT status
    if ~strcmp(solver.stats.return_status, 'Solve_Succeeded')
        % IPOPT fails
        disp('IPOPT fails')
        break
    else
        % pring information (solver.stats)
        msg = ['Homotopy Iter: ', num2str(homotopy_counter), '; ',...
            's: ', num2str(args.p,'%10.2e'), '; ',...
            'Cost: ', num2str(full(solution.f),'%10.2e'), '; ',...
            'Ipopt IterNum: ', num2str(solver.stats.iter_count), '; ',...
            'Time: ', num2str(1000*HomotopyTime,'%10.2e'), ' ms'];
        disp(msg)
    end
    % check termination of homotopy
    if (args.p == sEnd) && (strcmp(solver.stats.return_status, 'Solve_Succeeded'))
        % success
        break
    else
        % update x0 and s for next homotopy iteration
        args.x0 = solution.x;
        s_trail = min([kappa_s_times .* args.p, args.p.^kappa_s_exp]);
        args.p = max([s_trail, sEnd]);
        homotopy_counter = homotopy_counter + 1;        
    end
end
toc(totalTime_Start)

%%
XTAUPW_Opt = reshape(full(solution.x), (x_Dim + p_Dim + w_Dim), nStages);
x_Opt   = XTAUPW_Opt(1                 : x_Dim, :);
p_Opt   = XTAUPW_Opt(1 + x_Dim         : x_Dim + p_Dim, :);
w_Opt   = XTAUPW_Opt(1 + x_Dim + p_Dim : end, :);

timeAxis = 0 : timeStep : nStages * timeStep;

f_Fun_map = f_Fun.map(nStages);
f_value = f_Fun_map(x_Opt, p_Opt);
f_value = full(f_value);

figure(112)
subplot(3,2,1)
plot(timeAxis, [InitState, x_Opt], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time(s)')
title('system state')

subplot(3,2,2)
plot(timeAxis(2:end), f_value, 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(3,2,3)
plot(timeAxis(2:end), p_Opt, 'k', 'LineWidth', 1.2)
legend('y') 
xlabel('time(s)')
title('smoothing function')

subplot(3,2,4)
plot(timeAxis(2:end), x_Opt, 'g', 'LineWidth', 1.2)
legend('z := x') 
xlabel('time(s)')
title('switch function')

subplot(3,2,5)
plot(timeAxis(2:end), x_Opt, 'k',...
    timeAxis(2:end), p_Opt, 'b', 'LineWidth', 1.2)
legend('z(K)', 'y(p)') 
xlabel('time(s)')
title('checking BVI')
