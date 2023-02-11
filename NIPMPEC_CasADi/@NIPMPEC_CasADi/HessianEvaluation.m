function Hessian = HessianEvaluation(self, Var, Jac, s, mode, FRP)
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here

switch mode
    case 'Regular'
        %% regular iteration routine
        switch self.Option.HessianApproximation
            case 'Exact'
                LAG_hessian = self.FunObj.LAG_hessian(Var.x, Var.p, Var.w,...
                    Var.sigma, Var.eta, Var.gamma, s);
                Hessian = full(LAG_hessian);
            case 'CostFunction'
                L_hessian = self.FunObj.L_hessian(Var.x, Var.p, Var.w);
                Hessian = full(L_hessian);
            case 'GaussNewton'
                LAGx = Jac.Lx - Var.sigma' * Jac.Gx + Var.eta' * Jac.Cx - Var.gamma' * Jac.PHIx;
                LAGp = Jac.Lp - Var.sigma' * Jac.Gp + Var.eta' * Jac.Cp - Var.gamma' * Jac.PHIp;
                LAGw = Jac.Lw - Var.sigma' * Jac.Gw + Var.eta' * Jac.Cw - Var.gamma' * Jac.PHIw;
                Hessian = [LAGx, LAGp, LAGw]' * [LAGx, LAGp, LAGw];
            otherwise
                error('specified method to compute Hessian is not supported')
        end
               
    case 'FRP'
        %% FRP iteration routine
        FRP_L_hessian = self.FunObj.FRP_L_hessian(Var.x, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
        Hessian = full(FRP_L_hessian);
        
end

end

