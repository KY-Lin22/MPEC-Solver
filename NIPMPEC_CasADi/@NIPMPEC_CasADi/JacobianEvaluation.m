function Jac = JacobianEvaluation(self, Var, s, mode, FRP)
%JacobianEvaluation
%   return struct Jac with fileds including the following values
% Lx Lp Lw
% Gx Gp Gw
% Cx Cp Cw
% PHIx PHIp PHIw
%

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
%
Jac = struct('Lx', full(Lx), 'Lp', full(Lp), 'Lw', full(Lw),...
    'Gx', full(Gx), 'Gp', full(Gp), 'Gw', full(Gw),...
    'Cx', full(Cx), 'Cp', full(Cp), 'Cw', full(Cw),...
    'PHIx', full(PHIx), 'PHIp', full(PHIp), 'PHIw', full(PHIw));
end

