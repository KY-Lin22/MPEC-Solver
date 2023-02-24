function Fun = FunctionEvaluation(self, Var, s, z, mode, FRP)
%FunctionEvaluation
%   return struct Fun with fileds including the following values
% L G C PHI
% PSIg 
% PSIgSigma_diagVec 
% PSIgG_diagVec
% PSIphi
% PSIphiGamma_diagVec 
% PSIphiPHI_diagVec
% 

% cost function
switch mode
    case 'Regular'
        L = self.FunObj.L(Var.x, Var.p, Var.w);
    case 'FRP'
        L = self.FunObj.FRP_L(Var.x, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end

% constraint
G = self.FunObj.G(Var.x, Var.p, Var.w);
C = self.FunObj.C(Var.x, Var.p, Var.w);
PHI = self.FunObj.PHI(Var.x, Var.p, Var.w, s);

% FB function and Jacobian for G
PSIg = self.FunObj.FB_G(Var.sigma, full(G), z);
[PSIgSigma_diagVec, PSIgG_diagVec] = self.FunObj.FB_G_grad(Var.sigma, full(G), z);

% FB function and Jacobian for PHI
PSIphi = self.FunObj.FB_PHI(Var.gamma, full(PHI), z);
[PSIphiGamma_diagVec, PSIphiPHI_diagVec] = self.FunObj.FB_PHI_grad(Var.gamma, full(PHI), z);

%
Fun = struct('L', full(L),...
    'G', full(G), 'C', full(C), 'PHI', full(PHI),...
    'PSIg', full(PSIg),...
    'PSIgSigma_diagVec', full(PSIgSigma_diagVec),...
    'PSIgG_diagVec', full(PSIgG_diagVec),...
    'PSIphi', full(PSIphi),...
    'PSIphiGamma_diagVec', full(PSIphiGamma_diagVec),...
    'PSIphiPHI_diagVec', full(PSIphiPHI_diagVec));
end

