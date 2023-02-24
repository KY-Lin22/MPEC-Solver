function [dY, Info] = SearchDirection(self, KKT_Residual, KKT_Matrix)
%UNTITLED31 Summary of this function goes here
%   Detailed explanation goes here
TimeStart = tic;

% KKT residual and matrix
T = [KKT_Residual.G_Fsb;...
     KKT_Residual.C_Fsb;...
     KKT_Residual.PHI_Fsb;...
     KKT_Residual.LAGxT;...
     KKT_Residual.LAGpT;...
     KKT_Residual.LAGwT];
J = KKT_Matrix.J;

switch self.Option.linearSystemSolver
    case 'linsolve_Sym_dense'
        % solve linear system using linsolve(J is dense) with 'SYM' option
        opts.SYM = true;
        dY = linsolve(J, -T, opts);     
    case 'mldivide_dense'
        % solve linear system using mldivide(J is dense)
        dY = J\(-T);
    case 'mldivide_sparse'
        % solve linear system using mldivide(J is sparse)
        dY = J\(-T);
    case 'pinv'
        % solve linear system using pinv
        dY = pinv(J)*(-T);
    otherwise
        error('specified method is not supported')
end

TimeElapsed = toc(TimeStart);
Info.Time = TimeElapsed;
end

