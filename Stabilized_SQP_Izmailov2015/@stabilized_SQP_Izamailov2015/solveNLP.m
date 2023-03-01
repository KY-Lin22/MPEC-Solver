function [x_Opt, Info] = solveNLP(self, x_Init)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;
NLP = self.NLP;
Option = self.Option;
FunObj = self.FunObj;

% very first primal and dual variable
x_0 = x_Init;
lambda_0 = zeros(NLP.hDim, 1); % here I initilize the very first dual variables with value which satisfies their estimate
mu_0 = zeros(NLP.gDim, 1);     % bounds (paragraph 3 in page 424) so that their first estimate can choose these value.

% initialize estimate for the first dual variable (here 0 is the value of outer loop counter k)
bar_lambda_0 = lambda_0; % \in [Option.bar_lambdaMin, Option.bar_lambdaMax]
bar_mu_0 = mu_0; % \in [0, Option.bar_muMax]

% initialize estimate for the first primal variable (here 0 is the value of outer loop counter k and inner loop counter j)
hat_x_0 = x_0;

%% SQP iteration rountie
% k: outer iteration counter, increase after sSQP(step 2 in ref.) and Aug-L iteration(step 6 and 7 in ref.)
% j: inner iteration counter, increase after current Aug-L subproblem does not have acceptable stationary property(step 6 in ref.)
% in k iteration, x_k is the given iterate x_{k}, x_kNext is the new iterate x_{k+1}

% initialize given iterate used in outer loop iteration
hat_x_k = hat_x_0;
bar_lambda_k = bar_lambda_0;
bar_mu_k = bar_mu_0;
r_k = Option.r0;
epsilon_k = Option.epsilon0;
sigma_k = Option.sigma0;
iterateType_k = 'Init';
innerLoopFlag = 'Success';
disp('********************************************************************')
% outer loop iteration
for k = 0 : Option.maxOuterIterNum
    % STEP 1: compute natural residual rho of given iterate (k)
    rho_k = FunObj.rho(hat_x_k, bar_lambda_k, bar_mu_k);
    rho_k = full(rho_k);
    % STEP 2: check outer loop termination (end of paragraph 1 in page 424)
    disp(['Iter: ', num2str(k), '; ',...
        'iterateType: ', iterateType_k, '; ',...
        'rho: ', num2str(rho_k)])
    disp('--------------------------------------------------------------------')
    if rho_k < Option.tol_rho
        % solver finds the optimal solution
        exitFlag = true;
        terminalStatus = 'Success';
    elseif (k == Option.maxOuterIterNum) || (strcmp(innerLoopFlag, 'Fail'))
        % solver fails to find the optimal solution
        exitFlag = true;
        terminalStatus = 'Fail';
    else
        exitFlag = false;
    end
    % STEP 3: check exitFlag
    if exitFlag
        timeElapsed = toc(timeStart);
        % return given iterate as solution
        x_Opt = hat_x_k;
        % organize output information
        Info = struct('iterNum', k, 'terminalStatus', terminalStatus, 'time', timeElapsed, ...
            'lambda', bar_lambda_k, 'mu', bar_mu_k,...
            'rho', rho_k, 'iterateType', iterateType_k);
        % print information
        disp('********************************************************************')
        disp('Done!')
        disp(['iterNum: ', num2str(Info.iterNum)])
        disp(['terminalStatus: ', Info.terminalStatus])
        disp(['elapsedTime: ', num2str(Info.time), ' seconds'])
        disp(['rho: ', num2str(Info.rho)])
        disp(['solution type: ', Info.iterateType])
        break
    end
    % STEP 4: inner loop iteration (Algorithm 1)
    [hat_x_kNext, bar_lambda_kNext, bar_mu_kNext, r_kNext, epsilon_kNext, sigma_kNext, iterateType_kNext, innerLoopFlag] = ...
        self.innerLoopIteration(hat_x_k, bar_lambda_k, bar_mu_k, r_k, epsilon_k, sigma_k, iterateType_k);
    % STEP 5: prepare for next outer iteration
    hat_x_k = hat_x_kNext;
    bar_lambda_k = bar_lambda_kNext;
    bar_mu_k = bar_mu_kNext;
    r_k = r_kNext;
    epsilon_k = epsilon_kNext;
    sigma_k = sigma_kNext;
    iterateType_k = iterateType_kNext;
end
end

