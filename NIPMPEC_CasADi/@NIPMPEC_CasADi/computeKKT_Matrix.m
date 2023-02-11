function KKT_Matrix = computeKKT_Matrix(self, Fun, Jac, Hessian)
%UNTITLED30 Summary of this function goes here
%   Detailed explanation goes here

Dim = self.MPEC.Dim;
% regular parameter
nu_G = self.Option.RegularParam.nu_G;
nu_J = self.Option.RegularParam.nu_J;
nu_H = self.Option.RegularParam.nu_H;

% Constraint Jacobian
Gvar = [Jac.Gx, Jac.Gp, Jac.Gw];
Cvar = [Jac.Cx, Jac.Cp, Jac.Cw];
PHIvar = [Jac.PHIx, Jac.PHIp, Jac.PHIw];

% D
d = -(diag(Fun.PSIgSigma) - nu_G * ones(Dim.sigma, 1))./(diag(Fun.PSIgG) - nu_G * ones(Dim.sigma, 1));
D = diag(d);

% E
e = -(diag(Fun.PSIphiGamma) - nu_G * ones(Dim.gamma, 1))./(diag(Fun.PSIphiPHI) - nu_G * ones(Dim.gamma, 1));
E = diag(e);

% assemble KKT matrix
KKT_Matrix.J = [D,                            zeros(Dim.sigma, Dim.eta),  zeros(Dim.sigma, Dim.gamma),  -Gvar;...
                zeros(Dim.eta, Dim.sigma),    -nu_J * eye(Dim.eta),       zeros(Dim.eta, Dim.gamma),    Cvar;...
                zeros(Dim.gamma, Dim.sigma),  zeros(Dim.gamma, Dim.eta),  E,                            -PHIvar;...
                -Gvar',                       Cvar',                      -PHIvar',                     Hessian + nu_H*eye(Dim.Z)]; 

end

