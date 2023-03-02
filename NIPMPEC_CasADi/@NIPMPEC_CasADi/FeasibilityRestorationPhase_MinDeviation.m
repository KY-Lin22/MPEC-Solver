function [Var_FRP, Info] = FeasibilityRestorationPhase_MinDeviation(self, Var_Ref, Fun_Ref, Jac_Ref, s, z)
%UNTITLED35 Summary of this function goes here
%   Detailed explanation goes here

TimeStart = tic;

Dim = self.MPEC.Dim;
Option = self.Option;
maxIterNum = Option.FRP.maxIterNum;
nu_Z = Option.FRP.nu_Z;
nu_M = Option.FRP.nu_M;
betaInit = Option.FRP.betaInit;

if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp(' ')
    disp('*---- Entering Feasibility Restoration Phase...----*')
end

% compute ZRef and ZWeight
ZRef = [Var_Ref.x; Var_Ref.p; Var_Ref.w];
ZWeight = nu_Z * min([ones(Dim.Z, 1), 1./abs(ZRef)], [], 2);
FRP = struct('ZRef', ZRef, 'ZWeight', ZWeight);

% init Var (reusing Var_Ref except the dual variables of equality-type constraint)
Var_Init = struct('x', Var_Ref.x, 'p', Var_Ref.p, 'w', Var_Ref.w,...
    'sigma', Var_Ref.sigma, 'eta', zeros(Dim.eta, 1), 'gamma', Var_Ref.gamma);

% init Fun (reusing Fun_Ref except the cost function)
Fun_Init = struct('L', 0,...% no deviation for the VarInit
    'G', Fun_Ref.G, 'C', Fun_Ref.C, 'PHI',Fun_Ref.PHI,...
    'PSIg', Fun_Ref.PSIg,  'PSIphi', Fun_Ref.PSIphi);

% init Jac
[FRP_Lx_Init, FRP_Lp_Init, FRP_Lw_Init] = self.FunObj.FRP_L_grad(Var_Ref.x, Var_Ref.p, Var_Ref.w, FRP.ZRef, FRP.ZWeight);
Jac_Init = struct('Lx', full(FRP_Lx_Init), 'Lp', full(FRP_Lp_Init), 'Lw', full(FRP_Lw_Init),...
    'Gx', Jac_Ref.Gx, 'Gp', Jac_Ref.Gp,...
    'Cx', Jac_Ref.Cx, 'Cp', Jac_Ref.Cp, 'Cw', Jac_Ref.Cw,...
    'PHIp', Jac_Ref.PHIp, 'PHIw', Jac_Ref.PHIw,...
    'PSIgSigma_diagVec', Jac_Ref.PSIgSigma_diagVec, 'PSIgG_diagVec', Jac_Ref.PSIgG_diagVec,...
    'PSIphiGamma_diagVec', Jac_Ref.PSIphiGamma_diagVec, 'PSIphiPHI_diagVec', Jac_Ref.PSIphiPHI_diagVec);

% init Feasibility
Feasibility_Init = norm([Fun_Init.PSIg; Fun_Init.C; Fun_Init.PSIphi], 1);

%% solving MPEC
% initialize FRP iteration routine (x: previous iterate; x_j: current iterate)
Var = Var_Init;
Fun = Fun_Init;
beta = betaInit;

FailureFlag.LineSearch = false;

% FRP iteration routine 
for j = 1 : maxIterNum + 1
    %% step 1: compute totalCost and Feasibility of previous iterate
    totalCost = Fun.L;
    Feasibility = norm([Fun.PSIg; Fun.C; Fun.PSIphi], 1);
    % print information of previous iterate
    if Option.printLevel == 2
        prevIterMsg = ['Iter: ', num2str(j - 1), '; ',...
            'Cost: ', num2str(totalCost, '%10.3e'), '; ',...
            'Feasibility: ', num2str(Feasibility, '%10.3e'), '; ',...
            'beta: ' num2str(beta, '%10.3e')];
        disp(prevIterMsg)
    end     
    
    %% step 2: Checking Termination
    terminalCond.Feasibility = (Feasibility <= nu_M * Feasibility_Init);    
    terminalCond.maxIterNum = (j == (maxIterNum + 1));
    terminalCond.LS = FailureFlag.LineSearch;      
    
    if terminalCond.Feasibility
        % FRP finds a less infeasibility solution
        exitFlag = true;
        FRP_FailureFlag = false;
    elseif terminalCond.maxIterNum || terminalCond.LS
        % FRP fails to find a less infeasibility solution
        exitFlag = true;
        FRP_FailureFlag = true;
    else
        exitFlag = false;        
    end    
    
    %% step 3: Checking Exit Flag
    if exitFlag
        % return solution
        if ~FRP_FailureFlag
            % return a less infeasibility solution (reusing Var except the dual variables of equality-type constraint)
            Var_FRP = struct('x', Var.x, 'p', Var.p, 'w', Var.w,...
                'sigma', Var.sigma, 'eta', zeros(Dim.eta, 1), 'gamma', Var.gamma);
            % create Fun_FRP (reusing Fun except cost function)
            L = self.FunObj.L(Var_FRP.x, Var_FRP.p, Var_FRP.w);
            Fun_FRP = struct('L', full(L),...
                'G', Fun.G, 'C', Fun.C, 'PHI',Fun.PHI,...
                'PSIg', Fun.PSIg, 'PSIphi', Fun.PSIphi);
            % termination message
            terminationMsg_1 = ['FRP returns a less infeasibility iterate after ', num2str(j - 1), ' iteration; '];
            terminationMsg_2 = ['Feasibility: ', num2str(Feasibility_Init), '(Init)-->', num2str(Feasibility), '(End); ',...
                'DeviationCost: ' num2str(0), '(Init)-->', num2str(totalCost), '(End); '];
        else
            % return Var_Ref and Fun_Ref
            Var_FRP = Var_Ref;
            Fun_FRP = Fun_Ref;
            % termination message
            terminationMsg_1 = 'FRP fails to return a less infeasibility iterate';
            terminationMsg_2 = ['failure condition are: ',...
                'maxIterNum', '(', mat2str(terminalCond.maxIterNum), '); ', 'LineSearch', '(', mat2str(terminalCond.LS), '):'];
        end
        % record and print information about solution
        FRP_TimeElapsed = toc(TimeStart);
        Info.Fun = Fun_FRP;
        Info.Time = FRP_TimeElapsed;
        Info.FailureFlag = FRP_FailureFlag;          
        if (Option.printLevel == 1) || (Option.printLevel == 2)
            disp(terminationMsg_1)
            disp(terminationMsg_2)
        end
        
        break       
    end
    
    %% step 4: Evaluate KKT Residual and Matrix of Previous Iterate
    % Jacobian and KKT residual
    if j == 1
        Jac = Jac_Init;
    else       
        Jac = self.JacobianEvaluation(Var, Fun, s, z, 'FRP', FRP);
    end
    [KKT_Residual, ~] = self.computeKKT_Residual_Error(Var, Fun, Jac);
    % Hessian and KKT matrix
    Hessian = self.HessianEvaluation(Var, Jac, s, 'FRP', FRP);
    KKT_Matrix = self.computeKKT_Matrix(Jac, Hessian);
    
    %% step 5: Search Direction Evaluation
    [dY_k, ~] = self.SearchDirection(KKT_Residual, KKT_Matrix);
    
    %% step 6: Merit Line Search 
    [Var_LS, Info_LS] = self.LineSearch_Merit(Var, Fun, Jac, beta, s, z,...
        [], [], dY_k, 'FRP', FRP);
    % check failure flag and determine new iterate and its function evaluation
    if Info_LS.FailureFlag
        FailureFlag.LineSearch = true;
        % use previous iterate       
        Var_k  = Var;
        Fun_k  = Fun;
        beta_k = beta;               
    else
        FailureFlag.LineSearch = false;
        % use MLS iterate     
        Var_k  = Var_LS;
        Fun_k  = Info_LS.Fun;
        beta_k = Info_LS.beta;                
    end    
    % prepare for next iteration
    Var = Var_k;
    beta = beta_k;
    Fun = Fun_k;      
end

if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('*--------------------------------------------------*')
end

end

