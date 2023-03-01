function Jac = JacobianEvaluation(self, Var, Fun, s, z, mode, FRP)
%JacobianEvaluation
%   return struct Jac with fileds including the following values
% Lx Lp Lw
% Gx Gp Gw
% Cx Cp Cw
% PHIx PHIp PHIw
% PSIgSigma_diagVec 
% PSIgG_diagVec
% PSIphiGamma_diagVec 
% PSIphiPHI_diagVec

% cost function Jacobian
switch mode
    case 'Regular'
        [Lx, Lp, Lw] = self.FunObj.L_grad(Var.x, Var.p, Var.w);
    case 'FRP'
        [Lx, Lp, Lw] = self.FunObj.FRP_L_grad(Var.x, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end

% constraint Jacobian
[Gx, Gp, Gw] = self.FunObj.G_grad(Var.x, Var.p, Var.w);
[Cx, Cp, Cw] = self.FunObj.C_grad(Var.x, Var.p, Var.w);
[PHIx, PHIp, PHIw] = self.FunObj.PHI_grad(Var.x, Var.p, Var.w, s);

% FB Jacobian for G and PHI
[PSIgSigma_diagVec, PSIgG_diagVec] = self.FunObj.FB_G_grad(Var.sigma, Fun.G, z);
[PSIphiGamma_diagVec, PSIphiPHI_diagVec] = self.FunObj.FB_PHI_grad(Var.gamma, Fun.PHI, z);

%
Jac = struct('Lx', full(Lx), 'Lp', full(Lp), 'Lw', full(Lw),...
    'Gx', full(Gx), 'Gp', full(Gp), 'Gw', full(Gw),...
    'Cx', full(Cx), 'Cp', full(Cp), 'Cw', full(Cw),...
    'PHIx', full(PHIx), 'PHIp', full(PHIp), 'PHIw', full(PHIw),...
    'PSIgSigma_diagVec', full(PSIgSigma_diagVec),...
    'PSIgG_diagVec', full(PSIgG_diagVec),...
    'PSIphiGamma_diagVec', full(PSIphiGamma_diagVec),...
    'PSIphiPHI_diagVec', full(PSIphiPHI_diagVec));
end

