function Hessian = HessianEvaluation(self, Var, Jac, s, mode, FRP)
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here
Option = self.Option;
FunObj = self.FunObj;

switch mode
    case 'Regular'
        %% regular iteration routine
        switch Option.HessianApproximation
            case 'Exact'
                LAG_hessian = FunObj.LAG_hessian(Var.x, Var.p, Var.w,...
                    Var.sigma, Var.eta, Var.gamma, s);
                Hessian = LAG_hessian;
            case 'CostFunction'
                L_hessian = FunObj.L_hessian(Var.x, Var.p, Var.w);
                Hessian = L_hessian;
            case 'GaussNewton'
                Hessian = [Jac.Lx, Jac.Lp, Jac.Lw]' * [Jac.Lx, Jac.Lp, Jac.Lw];
            otherwise
                error('specified method to compute Hessian is not supported')
        end
               
    case 'FRP'
        %% FRP iteration routine
        FRP_L_hessian = FunObj.FRP_L_hessian(Var.x, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
        Hessian = FRP_L_hessian;
        
end

if strcmp(Option.linearSystemSolver, 'mldivide_sparse')
    Hessian = sparse(Hessian);
else
    Hessian = full(Hessian);
end
end

