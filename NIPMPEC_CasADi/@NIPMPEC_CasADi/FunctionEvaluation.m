function Fun = FunctionEvaluation(self, Var, s, z, mode, FRP)
%FunctionEvaluation
%   return struct Fun with fileds including the following values
% L G C PHI
% PSIg 
% PSIphi
% 
FunObj = self.FunObj;
% cost function
switch mode
    case 'Regular'
        L = FunObj.L(Var.x, Var.p, Var.w);
    case 'FRP'
        L = FunObj.FRP_L(Var.x, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end

% constraint
G = FunObj.G(Var.x, Var.p);
C = FunObj.C(Var.x, Var.p, Var.w);
PHI = FunObj.PHI(Var.p, Var.w, s);

% FB function for G and PHI
PSIg = FunObj.FB_G(Var.sigma, full(G), z);
PSIphi = FunObj.FB_PHI(Var.gamma, full(PHI), z);

%
Fun = struct('L', full(L),...
    'G', full(G), 'C', full(C), 'PHI', full(PHI),...
    'PSIg', full(PSIg), 'PSIphi', full(PSIphi));
end

