clc
clear all
%% problem formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 400;
sInit = 1e-1;
sEnd = 1e-8; 
zEnd = 1e-8;

x_Dim = 4;
tau_Dim = 1;
p_Dim = 1;
w_Dim = p_Dim; % auxiliary var
s_Dim = 1; % regularization parameter
x = SX.sym('x', x_Dim);
tau = SX.sym('tau', tau_Dim);
p = SX.sym('p', p_Dim);
w = SX.sym('w', w_Dim);
s = SX.sym('s', s_Dim);

InitState = [1; 0/180*pi; 0; 0];
RefState = [1; 180/180*pi; 0; 0];
EndState = [1; 180/180*pi; 0; 0];

% dynamics
mass = [1; 0.1];
linkLength = 1;
g = 9.8;
M_Mat = [mass(1) + mass(2),                 mass(2) * linkLength * cos(x(2));...
         mass(2) * linkLength * cos(x(2)),  mass(2) * linkLength^2];
C_Mat = [0,   -mass(2) * linkLength * x(4) * sin(x(2));...
         0,   0]; 
G_Mat = [0;...
       -mass(2) * g * linkLength * sin(x(2))];
Bu = [tau(1);...
      0]; 
P = [p(1);...
     0]; % friction bewteen cart and ground  
H = G_Mat + Bu + P - C_Mat * [x(3); x(4)];   
f = [x(3:4);...
    inv(M_Mat)*H];
f_Fun = Function('f_Fun',{x, tau, p}, {f}, {'x', 'tau', 'p'}, {'f'});

% equilibrium constraint
eqlbm.l = -2;
eqlbm.u = 2;
eqlbm.K = x(3);
K_Fun = Function('K_Fun', {x, tau, p}, {eqlbm.K}, {'x', 'tau', 'p'}, {'K'});

% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
BVI = [p - eqlbm.l;...
    eqlbm.u - p;...%  p in [l, u]   
    w - eqlbm.K;...% auxiliary variable
    s - (p - eqlbm.l) * w;...
    s + (eqlbm.u - p) * w];% regularization
BVI_Fun = Function('BVI_Fun', {x, tau, p, w, s}, {BVI}, {'x', 'tau', 'p', 'w', 's'}, {'BVI'});
lbg_BVI = [0; 0; 0; 0; 0];
ubg_BVI = [Inf; Inf; 0; Inf; Inf];

% cost function
xWeight_stage = [1; 100; 1; 1];
tauWeight_stage = 1 * ones(tau_Dim, 1);
pWeight_stage = 0.001 * ones(p_Dim, 1);
wWeight_stage = 0.001 * ones(w_Dim, 1);

L_stageCost = 0.5*(x - RefState)'*diag(xWeight_stage)*(x - RefState)...
    + 0.5*tau'*diag(tauWeight_stage)*tau...
    + 0.5*p'*diag(pWeight_stage)*p...
    + 0.5*w'*diag(wWeight_stage)*w;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, tau, p, w}, {L_stageCost},...
    {'x', 'tau', 'p', 'w'}, {'L_stageCost'});

xWeight_terminal = [1; 100; 10; 20];
L_terminalCost = 0.5*(x - EndState)'*diag(xWeight_terminal)*(x - EndState);
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, p}, {L_terminalCost}, {'x', 'p'}, {'L_terminalCost'});

% formulate NLP 
X = SX.sym('X', x_Dim, nStages); % state variable
TAU = SX.sym('TAU', tau_Dim, nStages); % control variable
P = SX.sym('P', p_Dim, nStages); % algebraic variable
W = SX.sym('W', w_Dim, nStages); % auxilary variable
S = SX.sym('S', s_Dim, 1);% regularizaton parameter

x_Max = [5; 240/180*pi; 20; 20];
x_Min = [0; -240/180*pi; -20; -20];
tau_Max = 30;
tau_Min = -30;

lbx = -Inf * ones(x_Dim + tau_Dim + p_Dim + w_Dim, nStages);
ubx = Inf * ones(x_Dim + tau_Dim + p_Dim + w_Dim, nStages);
lbx(1 : x_Dim + tau_Dim, :) = repmat([x_Min; tau_Min], 1, nStages);
ubx(1 : x_Dim + tau_Dim, :) = repmat([x_Max; tau_Max], 1, nStages);

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
    tau_n = TAU(:, n);
    p_n = P(:, n);
    w_n = W(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, tau_n, p_n, w_n);
    L = L + L_n*timeStep;    
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, p_n);
        L = L + L_terminal;
    end 
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, tau_n, p_n);   
    % reformulated equilibrium constraint
    BVI_n = BVI_Fun(x_n, tau_n, p_n, w_n, S); 
    % constraint function
    g(:, n) = [F_n;...
        BVI_n];     
end

OCPEC = struct('x', reshape([X; TAU; P; W], (x_Dim + tau_Dim + p_Dim + w_Dim) * nStages, 1),...
    'f', L,...    
    'g', reshape(g, g_Dim * nStages, 1),...
    'p', S);

args.lbx = reshape(lbx, [], 1);
args.ubx = reshape(ubx, [], 1);
args.lbg = reshape(lbg, [], 1);
args.ubg = reshape(ubg, [], 1);

% initialize x0 and s
args.x0 = zeros((x_Dim + tau_Dim + p_Dim + w_Dim) * nStages, 1);
args.p = sInit;

% option (solver.print_options)
Option = struct;
Option.print_time = false;
Option.ipopt.max_iter = 2000;
Option.ipopt.tol = 1e-6;
Option.ipopt.mu_target = 0.5 * (zEnd)^2;
Option.ipopt.print_level = 0;

% create solver
solver = nlpsol('solver', 'ipopt', OCPEC, Option);

%% homotopy
kappa_s_times = 0.2;
kappa_s_exp = 1.5;
homotopy_counter = 1;
totalIter = 0;
totalTime_Start = tic;
while true
    % Homotopy (outer) iteration
    HomotopyTime_Start = tic;
    solution = solver('x0', args.x0, 'p',args.p,...
        'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg);
    HomotopyTime = toc(HomotopyTime_Start);
    % check IPOPT status
    if strcmp(solver.stats.return_status, 'Solve_Succeeded')
        % pring information (solver.stats)
        totalIter = totalIter + solver.stats.iter_count;
        msg = ['Homotopy Iter: ', num2str(homotopy_counter), '; ',...
            's: ', num2str(args.p,'%10.2e'), '; ',...
            'Cost: ', num2str(full(solution.f),'%10.2e'), '; ',...
            'Ipopt IterNum: ', num2str(solver.stats.iter_count), '; ',...
            'Time: ', num2str(1000*HomotopyTime,'%10.2e'), ' ms'];
        disp(msg)
    else
        % IPOPT fails
        disp('IPOPT fails')
        break
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
disp(['total iteration: ', num2str(totalIter)])

%%
XTAUPW_Opt = reshape(full(solution.x), (x_Dim + tau_Dim + p_Dim + w_Dim), nStages);
x_Opt   = XTAUPW_Opt(1                           : x_Dim, :);
tau_Opt = XTAUPW_Opt(1 + x_Dim                   : x_Dim + tau_Dim, :);
p_Opt   = XTAUPW_Opt(1 + x_Dim + tau_Dim         : x_Dim + tau_Dim + p_Dim, :);
w_Opt   = XTAUPW_Opt(1 + x_Dim + tau_Dim + p_Dim : end, :);

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