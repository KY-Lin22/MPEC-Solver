function [solution, Info] = solveMPEC(self, Var_Init)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

%% check input

%%
Option = self.Option;
maxIterNum = Option.maxIterNum;
Tol = Option.Tolerance;
betaInit = Option.LineSearch.betaInit;
employFRP = Option.employFeasibilityRestorationPhase;

sInit = Option.sInit;
zInit = Option.zInit;
sEnd = Option.sEnd;
zEnd = Option.zEnd;

% create Record
Record = struct('s', zeros(maxIterNum + 1, 1), 'z', zeros(maxIterNum + 1, 1), 'VarType', [],...
    'merit', zeros(maxIterNum + 1, 1), 'beta', zeros(maxIterNum + 1, 1), 'stepSize', zeros(maxIterNum + 1, 1), ...
    'totalCost', zeros(maxIterNum + 1, 1), 'KKT_Error', [], 'Time', [],...
    'iterNum', 0, 'terminalStatus',0, 'terminalCond', []); 
Record.VarType = cell(maxIterNum + 1, 1);
Record.KKT_Error = struct('Total', zeros(maxIterNum + 1, 1), ...
    'Feasibility', zeros(maxIterNum + 1, 1), 'Stationarity', zeros(maxIterNum + 1, 1));   
Record.Time = struct('JacobianHessian', 0, 'KKT', 0, 'SearchDirection', 0, 'LineSearch', 0, 'FRP', 0,...
    'else', 0, 'total', 0);

%% solving MPEC
disp('******************************************************************');
disp('Computing the Optimal Solution for MPEC...');

% initialize regular iteration routine (x: previous iterate; x_k: current iterate)
VarType = 'Init';
Var = struct('x', Var_Init.x, 'p', Var_Init.p, 'w', Var_Init.w,...
    'sigma', Var_Init.sigma, 'eta', Var_Init.eta, 'gamma', Var_Init.gamma);
s = sInit;
z = zInit;
Fun = self.FunctionEvaluation(Var, s, z, 'Regular', []); 
merit = 0;
beta = betaInit;
stepSize = 0;

TimeElasped_JacobianHessian = 0;
TimeElasped_KKT = 0;
TimeElasped_SearchDirection = 0;
TimeElasped_LineSearch = 0;
TimeElasped_FRP = 0;
TimeElasped_Total = 0;

FailureFlag.LineSearch = false;
FailureFlag.FRP = false;

% regular iteration routine
for k = 1 : maxIterNum + 1
    totalTimeStart = tic;
    %% step 1: Evaluate KKT Residual of Previous Iterate
    % Jacobian
    Jacobian_TimeStart = tic;
    Jac = self.JacobianEvaluation(Var, s, 'Regular', []);
    TimeElasped_Jacobian = toc(Jacobian_TimeStart);
    
    % totalCost, KKT residual and error
    KKT_Residual_TimeStart = tic;
    totalCost = Fun.L;
    [KKT_Residual, KKT_Error] = self.computeKKT_Residual_Error(Var, Fun, Jac);
    TimeElasped_KKT_Residual = toc(KKT_Residual_TimeStart);
    
    %% Record and Print Information of Previous Iterate
    Record.s(k) = s;
    Record.z(k) = z;
    Record.VarType{k} = VarType;    
    
    Record.merit(k) = merit;
    Record.beta(k) = beta;
    Record.stepSize(k) = stepSize;
    
    Record.totalCost(k)              = totalCost;
    Record.KKT_Error.Total(k)        = KKT_Error.Total;
    Record.KKT_Error.Feasibility(k)  = KKT_Error.Feasibility;
    Record.KKT_Error.Stationarity(k) = KKT_Error.Stationarity;      

    TimeElasped_Else = TimeElasped_Total - TimeElasped_JacobianHessian - TimeElasped_KKT...
        - TimeElasped_SearchDirection - TimeElasped_LineSearch - TimeElasped_FRP;    
    Record.Time.JacobianHessian = Record.Time.JacobianHessian + TimeElasped_JacobianHessian;
    Record.Time.KKT             = Record.Time.KKT             + TimeElasped_KKT;
    Record.Time.SearchDirection = Record.Time.SearchDirection + TimeElasped_SearchDirection;
    Record.Time.LineSearch      = Record.Time.LineSearch      + TimeElasped_LineSearch;
    Record.Time.FRP             = Record.Time.FRP             + TimeElasped_FRP;    
    Record.Time.else            = Record.Time.else            + TimeElasped_Else;
    Record.Time.total           = Record.Time.total           + TimeElasped_Total;
    
    if Option.printLevel == 2
        prevIterMsg = ['Iter: ', num2str(k - 1), '(', VarType, ')', '; ',...
            's: ' num2str(s,'%10.2e'), '; ',...
            'z: ', num2str(z,'%10.2e'), '; ',...
            'Cost: ', num2str(totalCost,'%10.2e'), '; ',...
            'KKT: ',num2str(KKT_Error.Feasibility,'%10.2e'), '(F) ',...
                    num2str(KKT_Error.Stationarity,'%10.2e'),'(S); ',...
            'merit: ', num2str(merit,'%10.2e'), '; ',...
            'beta: ' num2str(beta,'%10.2e'), '; ',...
            'StepSize: ', num2str(stepSize,'%10.2e'), '; ', ...
            'Time: ', num2str(1000 * TimeElasped_Total,'%10.2e'), ' ms'];
        disp(prevIterMsg);
    end
    
    %% step 2: Checking Termination
    terminalCond.sz = ((s == sEnd) && (z == zEnd));
    terminalCond.KKT_T = (KKT_Error.Total <= Tol.KKT_Error_Total);
    terminalCond.KKT_F = (KKT_Error.Feasibility <= Tol.KKT_Error_Feasibility); % later transfer to (F <= Tol) && (active iteration counter)
    terminalCond.KKT_S = (KKT_Error.Stationarity <= Tol.KKT_Error_Stationarity);     
    terminalCond.maxIterNum = (k == (maxIterNum + 1));
    terminalCond.LSwithoutFRP = (FailureFlag.LineSearch && ~employFRP);
    terminalCond.LSwithFRP = (FailureFlag.LineSearch && FailureFlag.FRP);    
    if  terminalCond.sz && (terminalCond.KKT_T || terminalCond.KKT_F || terminalCond.KKT_S)
        % solver finds the optimal solution
        exitFlag = true;
        terminalStatus = 1;
    elseif terminalCond.maxIterNum || terminalCond.LSwithoutFRP || terminalCond.LSwithFRP
        % solver fails to find the optimal solution
        exitFlag = true;
        terminalStatus = 0;
    else
        exitFlag = false;
    end
    
    %% step 3: Checking Exit Flag
    if exitFlag
        % return solution (i.e.,iterate)
        Record.iterNum = k - 1; 
        Record.terminalStatus = terminalStatus;     
        Record.terminalCond = terminalCond;
        solution = Var; 
        Info = self.solutionExaminer(solution, Record);
        break
    end    
    
    %% step 4: Evaluate KKT Matrix of Previous Iterate
    % Hessian
    Hessian_TimeStart = tic;   
    Hessian = self.HessianEvaluation(Var, Jac, s, 'Regular', []);    
    TimeElasped_Hessian = toc(Hessian_TimeStart);
    TimeElasped_JacobianHessian = TimeElasped_Jacobian + TimeElasped_Hessian;
    % KKT matrix
    KKT_Matrix_TimeStart = tic;
    KKT_Matrix = self.computeKKT_Matrix(Fun, Jac, Hessian);
    TimeElasped_KKT_Matrix = toc(KKT_Matrix_TimeStart);
    TimeElasped_KKT = TimeElasped_KKT_Residual + TimeElasped_KKT_Matrix;
    
    %% step 5: Search Direction Evaluation
    [dY_k, Info_SD] = self.SearchDirection(KKT_Residual, KKT_Matrix);
    TimeElasped_SearchDirection = Info_SD.Time;
    
    %% step 6: Merit Line Search    
    [Var_LS, Info_LS] = self.LineSearch_Merit(Var, Fun, Jac, beta, s, z,...
        KKT_Residual, KKT_Matrix, dY_k, 'Regular', []);
    TimeElasped_LineSearch = Info_LS.Time;
    
    % check failure flag
    if Info_LS.FailureFlag
        FailureFlag.LineSearch = true;
        VarType_k = 'Prev';       
    else
        FailureFlag.LineSearch = false;
        VarType_k = Info_LS.VarType;      
    end 
    
    %% Step 7: Feasibility Restoration Phase
    if FailureFlag.LineSearch && employFRP
        [Var_FRP, Info_FRP] = self.FeasibilityRestorationPhase_MinDeviation(Var, Fun, Jac, s, z);
        TimeElasped_FRP = Info_FRP.Time;
        % check failure flag
        if Info_FRP.FailureFlag
            FailureFlag.FRP = true;
            VarType_k = 'Prev';
        else
            FailureFlag.FRP = false;           
            FailureFlag.LineSearch = false; % reset FailureFlag due to the success of FRP
            VarType_k = 'FRP';
        end        
    else       
        TimeElasped_FRP = 0;
    end       
    
    %% Step 8: Determine New Iterate and Its Function Evaluation
    if (strcmp(VarType_k,'MLS')) || (strcmp(VarType_k, 'SOC'))
        % use MLS or SOC Var
        Var_k  = Var_LS;   
        Fun_k  = Info_LS.Fun;
        merit_k    = Info_LS.merit;
        beta_k     = Info_LS.beta;
        stepSize_k = Info_LS.stepSize;               
    elseif (strcmp(VarType_k,'FRP'))
        % use FRP Var
        Var_k  = Var_FRP;  
        Fun_k  = Info_FRP.Fun;
        merit_k    = 0;
        beta_k     = betaInit;% reset penalty parameter
        stepSize_k = 1;                  
    elseif (strcmp(VarType_k,'Prev'))
        % use previous Var
        Var_k  = Var;
        Fun_k  = Fun;
        merit_k    = merit;
        beta_k     = beta;
        stepSize_k = stepSize;       
    else
        error('wrong VarType')
    end           
    
    %% Step 9: Update Perturbed Parameter and Prepare for Next Iteration
    % update perturbed parameter s, z and Fun_k
    if (strcmp(VarType_k,'MLS')) || (strcmp(VarType_k, 'SOC')) || (strcmp(VarType_k,'FRP'))
        [s_k, z_k, Fun_k] = self.computePerturedParam(Var_k, Fun_k, s, z);
    elseif  (strcmp(VarType_k,'Prev'))
        s_k = s;
        z_k = z;        
    else
        error('wrong VarType')
    end    
    % reset penalty parameter
    if (s_k ~= s) || (z_k ~= z)        
        beta_k = betaInit;
    end      
    % prepare for next iteration
    VarType = VarType_k;
    Var = Var_k;        
    s = s_k;
    z = z_k;  
    Fun = Fun_k;
    merit = merit_k;
    beta = beta_k;
    stepSize = stepSize_k;   
    
    TimeElasped_Total = toc(totalTimeStart);
  
end

disp('******************************************************************');

end

