function Option = createOption(self)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Option = struct;

Option.r0 = 1e4; % r0 > 0, init value of natural residual tolerance, used and updated in step 2, 7,
% default: 1e4, paragraph 2 in page 424
Option.epsilon0 = 1e2; % epsilon0 > 0, init value of augmented Lagrangian stationary tolerance, used in step 6, updated in step 2,7,
% default: 1e2, paragraph 2 in page 424
Option.sigma0 = 1e-4; % sigma0 > 0, init value of dual stabilization parameter, used in step 1(QP), 3,5,6,7(AugL), updated in step 2,7
% default: 1e-4, paragraph 3 in page 424
Option.gamma = 1; % gamma > 0, weight parameter in equ(11) checking the descent for AugL, used in the step 3
% default: 1, paragraph 3 in page 424

Option.q = 0.5; % 0 < q < 1, weight parameter in equ(10) updating natural residual tolerance rk, used in step 2, 7
% default: 0.5, paragraph 2 in page 424
Option.xita = 0.5; % 0 < xita < 1, weight parameter in step 7 updating augmented Lagrangian stationary tolerance epsilonk, used in step 7
% default: 0.5, paragraph 2 in page 424
Option.tau = 0.5; % 0 < tau < 1, weight parameter for step size in line search, used in step 5
% default: 0.5,  paragraph 3 in page 424
Option.epsilon = 0.1; % 0 < epsilon < 1, weight parameter for merit function in line search, used in step 5
% default: 0.1,  paragraph 3 in page 424
Option.kappa = 0.1; % 0 < kappa < 1, weight parameter in equ(15) for updating dual stabilization parameter sigma_k, used in step 7
% default: 0.1,  paragraph 3 in page 424
Option.delta = 0.5; % 0 < delta < 1, weight parameter in equ(15) for updating dual stabilization parameter sigma_k, used in step 7
% default: 0.5,  paragraph 3 in page 424

% bounds to safeguard the dual iterate
Option.bar_lambdaMin = -1e10; % paragraph 3 in page 424
Option.bar_lambdaMax = 1e10; % paragraph 3 in page 424
Option.bar_muMax = 1e10; % paragraph 3 in page 424

% modified Hessian
Option.omegaMax = 1e10;

% termination condition
Option.maxOuterIterNum = 500; % kMax, positive int(default: 500, end of paragraph 1 in page 424)
Option.maxInnerIterNum = 50; % jMax, positive int(default: none, so I specify by myself, although
% proposition 1 shows that j must be finite)
Option.maxLineSearchIterNum = 8; %iMax, nonnegative int(default: none, so I specify by myself)
Option.tol_rho = 1e-6; % tolerance for natural residual (default: 1e-6, end of paragraph 1 in page 424)
end

