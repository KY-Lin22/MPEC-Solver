function [Var_LS, Info] = LineSearch_Merit(self, Var, Fun, Jac, beta, s, z,...
    KKT_Residual, KKT_Matrix, dY_k, mode, FRP)
%UNTITLED32 Summary of this function goes here
%   Detailed explanation goes here

TimeStart = tic;

Dim = self.MPEC.Dim;
Option = self.Option;
employSOC = Option.employSecondOrderCorrection;

nu_D = Option.LineSearch.nu_D;
stepSize_Init = 1;
switch mode
    case 'Regular'
        rho = Option.LineSearch.rho;
        stepSize_Min = Option.LineSearch.stepSize_Min;
        stepSize_DecayRate = Option.LineSearch.stepSize_DecayRate;
    case 'FRP'
        rho = Option.FRP.rho;
        stepSize_Min = Option.FRP.stepSize_Min;
        stepSize_DecayRate = Option.FRP.stepSize_DecayRate;      
end

%% Some Evaluation Quantities of Previous Iterate
% cost and its directional derivative
totalCost = Fun.L;
LZ = [Jac.Lx, Jac.Lp, Jac.Lw];
dZ_k = dY_k(Dim.Node(3) + 1 : Dim.Node(6), :);
totalCostDD = LZ * dZ_k;

% constraint violation (L1 norm)
totalCstrVio_L1Norm = norm([Fun.PSIg; Fun.C; Fun.PSIphi], 1);

% penalty parameter
beta_Trial = totalCostDD/((1 - rho) * totalCstrVio_L1Norm);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end

% merit and its directional derivative 
merit = totalCost + beta_k * totalCstrVio_L1Norm;
meritDD = totalCostDD - beta_k * totalCstrVio_L1Norm;

%% Backtracking Line Search
hasFoundNewIterate = false;

while ~hasFoundNewIterate
    %% Step 1: estimate trail stepsize, iterate and merit 
    stepSize_trial = max([stepSize_Init, stepSize_Min]);
    Var_trial.sigma  = Var.sigma  + stepSize_trial * dY_k(              1 : Dim.Node(1), :);
    Var_trial.eta    = Var.eta    + stepSize_trial * dY_k(Dim.Node(1) + 1 : Dim.Node(2), :);
    Var_trial.gamma  = Var.gamma  + stepSize_trial * dY_k(Dim.Node(2) + 1 : Dim.Node(3), :);
    Var_trial.x      = Var.x      + stepSize_trial * dY_k(Dim.Node(3) + 1 : Dim.Node(4), :);
    Var_trial.p      = Var.p      + stepSize_trial * dY_k(Dim.Node(4) + 1 : Dim.Node(5), :);
    Var_trial.w      = Var.w      + stepSize_trial * dY_k(Dim.Node(5) + 1 : Dim.Node(6), :);
    
    Fun_trial = self.FunctionEvaluation(Var_trial, s, z, mode, FRP);
    totalCost_trail = Fun_trial.L;
    totalCstrVio_L1Norm_trial = norm([Fun_trial.PSIg; Fun_trial.C; Fun_trial.PSIphi], 1);   
    merit_trial = totalCost_trail + beta_k * totalCstrVio_L1Norm_trial;
    
    %% Step 2: checking sufficient decrease condition
    if merit_trial <= merit + stepSize_trial * nu_D * meritDD
        % return merit line search Var
        hasFoundNewIterate = true;
        LineSearchFailureFlag = false;
        VarType = 'MLS';         
        
    elseif (strcmp(mode, 'Regular')) && (employSOC) && (stepSize_trial == 1) && (totalCost_trail <= totalCost)
        % estimate second order correction (SOC) Var and merit
        Fun_full = Fun_trial;
        Var_trial = self.SecondOrderCorrection(Var, Fun, KKT_Residual, KKT_Matrix, Fun_full);
        
        Fun_trial = self.FunctionEvaluation(Var_trial, s, z, 'Regular', []);
        totalCost_trail = Fun_trial.L;
        totalCstrVio_L1Norm_trial = norm([Fun_trial.PSIg; Fun_trial.C; Fun_trial.PSIphi], 1);
        merit_trial = totalCost_trail + beta_k * totalCstrVio_L1Norm_trial;
        
        if merit_trial <= merit + stepSize_trial * nu_D * meritDD
            % return SOC Var
            hasFoundNewIterate = true;
            LineSearchFailureFlag = false;
            VarType = 'SOC'; 
        else
            % discard this SOC step and resume backtracking line search with smaller stepsize
            stepSize_Init = stepSize_DecayRate * stepSize_Init;
        end
    else
        % need to estimate a smaller stepsize
        stepSize_Init = stepSize_DecayRate * stepSize_Init;
    end  
    
    %% Step 3: checking min stepsize
    if (stepSize_trial == stepSize_Min)&&(~hasFoundNewIterate)
        LineSearchFailureFlag = true;
        VarType = 'MLS'; % return Var obtained by the stepSize_Min   
        break
    end     
    
end

%% Record Information
TimeElapsed = toc(TimeStart);

Info.VarType = VarType;
Var_LS = Var_trial;
Info.Fun = Fun_trial;
Info.merit = merit_trial;
Info.beta = beta_k;
Info.stepSize = stepSize_trial;

Info.FailureFlag = LineSearchFailureFlag;
Info.Time = TimeElapsed;

end
