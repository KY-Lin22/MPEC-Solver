function generateInitialGuess(self)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

Dim = self.MPEC.Dim;
% init
Var = struct('x', [], 'p', [], 'w', [],...
    'sigma', [], 'eta', [], 'gamma', []);

disp('Generating Initial Guess...')

%% generate initial guess for primal variable
Var.x = zeros(Dim.x, 1);
Var.p = 1/2 * (self.MPEC.l + self.MPEC.u);
K = self.FunObj.K(Var.x, Var.p);
Var.w = full(K);

%% generate initial guess for dual variable
Var.sigma  = ones(Dim.sigma, 1); % sigma >= 0
Var.eta = randn(Dim.eta, 1);
Var.gamma = ones(Dim.gamma, 1); % gamma >= 0

%% save initial guess
save('Gen_InitialGuess.mat', 'Var');
disp('Done!')

end

