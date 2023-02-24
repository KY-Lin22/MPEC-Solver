function [KKT_Residual, KKT_Error] = computeKKT_Residual_Error(self, Var, Fun, Jac)
%UNTITLED27 Summary of this function goes here
%   Detailed explanation goes here

Dim = self.MPEC.Dim;
nu_G = self.Option.RegularParam.nu_G;

%% KKT Residual
% primal feasibility
G_Fsb = - Fun.PSIg./(Fun.PSIgG_diagVec - nu_G * ones(Dim.sigma, 1));
C_Fsb = Fun.C;
PHI_Fsb = - Fun.PSIphi./(Fun.PSIphiPHI_diagVec - nu_G * ones(Dim.gamma, 1));

% dual feasibility
LAGx = Jac.Lx - Var.sigma' * Jac.Gx + Var.eta' * Jac.Cx - Var.gamma' * Jac.PHIx;
LAGp = Jac.Lp - Var.sigma' * Jac.Gp + Var.eta' * Jac.Cp - Var.gamma' * Jac.PHIp;
LAGw = Jac.Lw - Var.sigma' * Jac.Gw + Var.eta' * Jac.Cw - Var.gamma' * Jac.PHIw;

%
KKT_Residual = struct('G_Fsb', G_Fsb, 'C_Fsb', C_Fsb, 'PHI_Fsb', PHI_Fsb,...
    'LAGxT', LAGx', 'LAGpT', LAGp', 'LAGwT', LAGw');

%% KKT Error
% compute scaling parameter for KKT stationarity
scaling_max = 100;
LAMBDA_norm = norm(Var.sigma, 1) + norm(Var.eta, 1) + norm(Var.gamma, 1);    
LAMBDA_trial = (LAMBDA_norm)/(Dim.LAMBDA);
stationarityScaling = max([scaling_max, LAMBDA_trial])/scaling_max;

% compute KKT_Error
Feasibility = [Fun.PSIg; Fun.C; Fun.PSIphi];
Stationarity = [KKT_Residual.LAGxT; KKT_Residual.LAGpT; KKT_Residual.LAGwT];

KKT_Error.Feasibility = norm(Feasibility, Inf);
KKT_Error.Stationarity = norm(Stationarity, Inf)/stationarityScaling;
KKT_Error.Total = max([KKT_Error.Feasibility, KKT_Error.Stationarity]);

end

