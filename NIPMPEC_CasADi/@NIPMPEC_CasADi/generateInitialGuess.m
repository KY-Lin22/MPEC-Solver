function generateInitialGuess(self)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
MPEC = self.MPEC;
FunObj = self.FunObj;
Dim = MPEC.Dim;
% init
Var = struct('x', [], 'p', [], 'w', [],...
    'sigma', [], 'eta', [], 'gamma', []);

disp('Generating Initial Guess...')

%% generate initial guess for primal variable
Var.x = zeros(Dim.x, 1);
p = zeros(Dim.p, 1);
for i = 1 : Dim.p
    if (MPEC.l(i) == 0) && (MPEC.u(i) == Inf)
        % nonlinear complementary problem
        p(i, 1) = 1; % p > 0
    else
        % box constraint variation inequality
        p(i, 1) = 1/2 * (MPEC.l(i) + MPEC.u(i)); % l < = p < = u
    end
end
Var.p = p;
K = FunObj.K(Var.x, Var.p);
Var.w = full(K);

%% generate initial guess for dual variable
Var.sigma  = ones(Dim.sigma, 1); % sigma >= 0
Var.eta = randn(Dim.eta, 1);
Var.gamma = ones(Dim.gamma, 1); % gamma >= 0

%% save initial guess
save('Gen_InitialGuess.mat', 'Var');
disp('Done!')

end

