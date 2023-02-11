function Var_SOC = SecondOrderCorrection(self, Var, Fun, KKT_Residual, KKT_Matrix, Fun_full)

Dim = self.MPEC.Dim;
nu_G = self.Option.RegularParam.nu_G;

% SOC KKT residual
KKT_Residual_SOC = struct('G_Fsb', [], 'C_Fsb', [], 'PHI_Fsb', [],...
    'LAGxT', KKT_Residual.LAGxT, 'LAGpT', KKT_Residual.LAGpT, 'LAGwT', KKT_Residual.LAGwT);
G_Fsb_Correct = - Fun_full.PSIg./(diag(Fun.PSIgG) - nu_G * ones(Dim.sigma, 1));
KKT_Residual_SOC.G_Fsb = KKT_Residual.G_Fsb + G_Fsb_Correct;

KKT_Residual_SOC.C_Fsb = KKT_Residual.C_Fsb + Fun_full.C;

PHI_Fsb_Correct = - Fun_full.PSIphi./(diag(Fun.PSIphiPHI) - nu_G * ones(Dim.gamma, 1));
KKT_Residual_SOC.PHI_Fsb = KKT_Residual.PHI_Fsb + PHI_Fsb_Correct;

% compute second order correction direction and var
[dY_SOC, ~] = self.SearchDirection(KKT_Residual_SOC, KKT_Matrix);

Var_SOC.sigma  = Var.sigma  + dY_SOC(              1 : Dim.Node(1), :);
Var_SOC.eta    = Var.eta    + dY_SOC(Dim.Node(1) + 1 : Dim.Node(2), :);
Var_SOC.gamma  = Var.gamma  + dY_SOC(Dim.Node(2) + 1 : Dim.Node(3), :);
Var_SOC.x      = Var.x      + dY_SOC(Dim.Node(3) + 1 : Dim.Node(4), :);
Var_SOC.p      = Var.p      + dY_SOC(Dim.Node(4) + 1 : Dim.Node(5), :);
Var_SOC.w      = Var.w      + dY_SOC(Dim.Node(5) + 1 : Dim.Node(6), :);
end

