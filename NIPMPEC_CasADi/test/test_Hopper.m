clear all
clc
delete Gen_InitialGuess.mat

%% problem formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 200;
x_Dim = 8;
tau_Dim = 3;
p_Dim = 3;

InitState = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
MidState1 = [0.3; 0.65; 0; 0.2; 0; 0; 0; 0];
MidState2 = [0.4; 0.5; 0; 0.5; 0; 0; 0; 0];
MidState3 = [0.6; 0.65; 0; 0.2; 0; 0; 0; 0];
EndState = [0.7; 0.5; 0; 0.5; 0; 0; 0; 0];

xRef_init_mid1 = TrajectoryInterpolation(InitState, MidState1, 40);
xRef_mid1_mid2 = TrajectoryInterpolation(MidState1, MidState2, 40);
xRef_mid2_mid3 = TrajectoryInterpolation(MidState2, MidState3, 40);
xRef_mid3_end = TrajectoryInterpolation(MidState3, EndState, 40);
xRef_end_end = TrajectoryInterpolation(EndState, EndState, 40);
RefState = [xRef_init_mid1, xRef_mid1_mid2, xRef_mid2_mid3, xRef_mid3_end, xRef_end_end];

x = SX.sym('x', x_Dim, 1);
tau = SX.sym('tau', tau_Dim, 1);
p = SX.sym('p', p_Dim, 1); 

% hopper configuration (kinematic) for animation
Hopper_X_formula = [x(1);...
    x(1) + x(4)*sin(x(3))];
Hopper_Y_formula = [x(2);...
    x(2) - x(4)*cos(x(3))];
HopperConfiguration_Fun = Function('HopperConfiguration_Fun', {x}, {Hopper_X_formula, Hopper_Y_formula}, {'x'}, {'Hopper_X', 'Hopper_Y'});

% dynamics
mb = 1;
ml = 0.1;
Ib = 0.25;
Il = 0.025;
mu = 1;
g = 9.8;

M = diag([mb + ml, mb + ml, Ib + Il, ml]);
C = [0;...
    (mb + ml)*g;...
    0;...
    0];% here C := C*dq - G
B = [0, -sin(x(3));...
    0, cos(x(3));...
    1, 0;...
    0, 1];
W_N = [0, 1, x(4)*sin(x(3)), -cos(x(3))];
W_T = [1, 0, x(4)*cos(x(3)), sin(x(3))];
H = -C + B * [tau(1); tau(2)] + W_N' * p(1) + W_T' * tau(3);
f = [x(5:8);...
    inv(M)*H];
f_Fun = Function('f_Fun', {x, tau, p}, {f}, {'x', 'tau', 'p'}, {'f'});

% equilibrium constraint
eqlbm.l = [0; 0; 0];
eqlbm.u = [Inf;Inf;Inf];
Gap = x(2) - x(4) * cos(x(3));
pT_lb = tau(3) - (-mu * p(1));
pT_ub = (mu * p(1)) - tau(3);
eqlbm.K = [Gap; pT_lb; pT_ub];
K_Fun = Function('K_Fun', {x, tau, p}, {eqlbm.K}, {'x', 'tau', 'p'}, {'K'});

% cost function
xWeight_stage = [100; 100; 100; 100; 0.1; 0.1; 0.1; 0.1];
tauWeight = [0.1; 0.1; 0.001];
xRef = SX.sym('xRef', x_Dim, 1);
L_stageCost = 0.5*(x - xRef)'*diag(xWeight_stage)*(x - xRef) + 0.5*tau'*diag(tauWeight)*tau;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, tau, p, xRef}, {L_stageCost}, {'x', 'tau', 'p', 'xRef'}, {'L_stageCost'});

xWeight_terminal = [100; 100; 100; 100; 0.1; 0.1; 0.1; 0.1];
L_terminalCost = 0.5*(x - EndState)'*diag(xWeight_terminal)*(x - EndState);
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, tau, p}, {L_terminalCost}, {'x', 'tau', 'p'}, {'L_terminalCost'});

% inequality constraint
vel_T = W_T * x(5:8);
x_Max = [0.8; 0.7; pi; 0.5; 10; 10; 5; 5];
x_Min = [0; 0; -pi; 0.2; -10; -10; -5; -5];
tau_Max = [50; 50; 100];
tau_Min = [-50; -50; -100];
G_formula =...
    [0.01 - p(1) * vel_T;...% penalty slip motion
    x_Max - x;...
    x - x_Min;...
    tau_Max - tau;...
    tau - tau_Min];
G_formula_Fun = Function('G_formula_Fun', {x, tau, p}, {G_formula}, {'x', 'tau', 'p'}, {'G'});

% equality constraint
C_formula = p(2) - p(3) - vel_T; % using two auxilary variable for vel_T to reformulate friction 
C_formula_Fun = Function('C_formula_Fun', {x, tau, p}, {C_formula}, {'x', 'tau', 'p'}, {'C'});

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
    xRef_n = RefState(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, tau_n, p_n, xRef_n);
    L = L + L_n*timeStep;
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, tau_n, p_n);
        L = L + L_terminal;
    end    
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, tau_n, p_n);
    % equality constraint
    C_n = C_formula_Fun(x_n, tau_n, p_n);
    C(:, n) = [F_n;...
        C_n];    
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
Var_Init.x = reshape([RefState; randn(tau_Dim, nStages)], [], 1);

%% solving MPEC
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-2;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-4;

solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;
solver.Option.linearSystemSolver = 'mldivide_sparse'; % 'linsolve_Sym_dense', 'mldivide_dense', 'mldivide_sparse', 'pinv'

solver.Option.LineSearch.stepSize_Min = 0.01;
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

XTAU_Opt = reshape(solution.x, (x_Dim + tau_Dim), nStages);
x_Opt = XTAU_Opt(1 : x_Dim, :);
tau_Opt = XTAU_Opt(x_Dim + 1 : end, :);
p_Opt = reshape(solution.p, p_Dim, nStages);
timeAxis = 0 : timeStep : nStages * timeStep;

K_Fun_map = K_Fun.map(nStages);
K_value = K_Fun_map(x_Opt, tau_Opt, p_Opt);
K_value = full(K_value);

figure(111)
subplot(3,1,1)
plot(timeAxis, [InitState(1), x_Opt(1,:)], 'r',...
     timeAxis, [InitState(2), x_Opt(2,:)], 'g',...
     timeAxis, [InitState(3), x_Opt(3,:)], 'b',...
     timeAxis, [InitState(4), x_Opt(4,:)], 'k', 'LineWidth',1.2)
legend('x_b', 'y_b', '\theta_b', 'l_l')
legend('Orientation','horizontal')
xlabel('time [s]')
title('position')

subplot(3,1,2)
plot(timeAxis, [InitState(5), x_Opt(5,:)], 'r',...
     timeAxis, [InitState(6), x_Opt(6,:)], 'g',...
     timeAxis, [InitState(7), x_Opt(7,:)], 'b',...
     timeAxis, [InitState(8), x_Opt(8,:)], 'k', 'LineWidth',1.2)
xlabel('time [s]')
title('velocity')

subplot(3,1,3)
plot(timeAxis(2:end), tau_Opt(1,:),'b',...
    timeAxis(2:end), tau_Opt(2, :), 'k','LineWidth', 1.2)
legend('\tau_b', '\tau_l')
legend('Orientation','horizontal')
xlabel('time [s]')
title('control')

figure(112)
subplot(4,1,1)
plot(timeAxis(2:end), K_value(1, :), 'b','LineWidth', 1.2)
xlabel('time [s]')
ylabel('gap [m]')
title('normal impact contact')
subplot(4,1,2)
plot(timeAxis(2:end), p_Opt(1, :), 'k','LineWidth', 1.2)
xlabel('time [s]')
ylabel('impact [N]')

subplot(4,1,3)
plot(timeAxis(2:end), (p_Opt(2, :) - p_Opt(3, :)), 'b', 'LineWidth', 1.2)
xlabel('time [s]')
ylabel('foot vel [m/s]')
title('tangential friction contact')
subplot(4,1,4)
plot(timeAxis(2:end), tau_Opt(3, :), 'k', 'LineWidth', 1.2)
xlabel('time [s]')
ylabel('friction[N]')

%% animation
% get base line
baseline_X = [x_Min(1, 1); x_Max(1, 1)];  
baseline_Y = [0;0];

% get hopper position sequence based on given x sequence
HopperConfiguration_Fun_map = HopperConfiguration_Fun.map(nStages + 1);
[Hopper_X, Hopper_Y] = HopperConfiguration_Fun_map([InitState, x_Opt]);
Hopper_X = full(Hopper_X);
Hopper_Y = full(Hopper_Y);
trajbody_X = cell(1, nStages + 1);
trajbody_Y = cell(1, nStages + 1);
trajfoot_X = cell(1, nStages + 1);
trajfoot_Y = cell(1, nStages + 1);
for n = 1 : nStages + 1
    timeAxisSequence = [timeAxis(1:n), repmat(timeAxis(n), 1, nStages + 1 - n)]; 
    trajbody_X{1, n} = [Hopper_X(1, 1:n), repmat(Hopper_X(1, n), 1, nStages + 1 -n)];
    trajbody_Y{1, n} = [Hopper_Y(1, 1:n), repmat(Hopper_Y(1, n), 1, nStages + 1 -n)];    
    trajfoot_X{1, n} = [Hopper_X(2, 1:n), repmat(Hopper_X(2, n), 1, nStages + 1 -n)];
    trajfoot_Y{1, n} = [Hopper_Y(2, 1:n), repmat(Hopper_Y(2, n), 1, nStages + 1 -n)];    
end

% get figure size
figure(100)
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
% define movie record
mov = zeros(height*3/2, width*3/2, 1, nStages + 1, 'uint8');
% pre allocate
subplot(1, 1, 1)
plot(baseline_X, baseline_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
hold on
% init, some middle and end state
MS_ID = [1; nStages + 1];
for i = 1 : size(MS_ID, 1)
    MSi = MS_ID(i);
    plot(Hopper_X(:, MSi), Hopper_Y(:, MSi), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1)
    hold on
end
% trajectory
Hopper = plot(Hopper_X(:, 1), Hopper_Y(:, 1), '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 15, 'LineWidth', 1);
hold on
trajbody = plot(trajbody_X{1, 1}, trajbody_Y{1, 1}, '.--g', 'MarkerSize', 1, 'LineWidth', 1);
hold on
trajfoot = plot(trajfoot_X{1, 1}, trajfoot_Y{1, 1}, '.--r', 'MarkerSize', 1, 'LineWidth', 1);
hold on
% axis equal
axisLimit_X = baseline_X;
axisLimit_Y = [-0.1; 0.7];
axis([axisLimit_X; axisLimit_Y]);
xlabel('x_b [m]')
ylabel('y_b [m]')
timePrint = title(sprintf('Time: %0.2f sec', timeAxis(1)));
% animate
for n = 1 : nStages + 1
    % update XData and YData
    set(timePrint, 'String', sprintf('Time: %0.2f sec', timeAxis(n)));
    set(Hopper, 'XData', Hopper_X(:, n), 'YData', Hopper_Y(:, n));
    set(trajbody, 'XData', trajbody_X{1, n}, 'YData', trajbody_Y{1, n});
    set(trajfoot, 'XData', trajfoot_X{1, n}, 'YData', trajfoot_Y{1, n});   
    % get frame as an image
    f = getframe(gcf);
    % Create a colormap for the first frame. for the rest of the frames, use the same colormap
    if n == 1
        [mov(:,:,1,n), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,n) = rgb2ind(f.cdata, map, 'nodither');
    end        
end
% create an animated GIF
imwrite(mov, map, 'Hopper.gif', 'DelayTime', 0, 'LoopCount', inf)