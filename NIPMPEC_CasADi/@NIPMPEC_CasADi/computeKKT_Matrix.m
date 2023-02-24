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

% diag matrix vector: D = diag(d), nu_J,  E = diag(e);
d = -(Fun.PSIgSigma_diagVec - nu_G * ones(Dim.sigma, 1))./(Fun.PSIgG_diagVec - nu_G * ones(Dim.sigma, 1));
nu_J_vec = -nu_J*ones(Dim.eta, 1);
e = -(Fun.PSIphiGamma_diagVec - nu_G * ones(Dim.gamma, 1))./(Fun.PSIphiPHI_diagVec - nu_G * ones(Dim.gamma, 1));
diagVec = [d; nu_J_vec; e];

%% assemble KKT matrix
%J = [D,                            zeros(Dim.sigma, Dim.eta),  zeros(Dim.sigma, Dim.gamma),  -Gvar;...
%     zeros(Dim.eta, Dim.sigma),    -nu_J * eye(Dim.eta),       zeros(Dim.eta, Dim.gamma),    Cvar;...
%     zeros(Dim.gamma, Dim.sigma),  zeros(Dim.gamma, Dim.eta),  E,                            -PHIvar;...
%     -Gvar',                       Cvar',                      -PHIvar',                     Hessian + nu_H*eye(Dim.Z)];      

J = zeros(Dim.Y, Dim.Y);
for i = 1 : Dim.Node(3)
    J(i, i) = diagVec(i);
end
J(1 : Dim.Node(3), Dim.Node(3) + 1 : Dim.Node(6)) = [-Gvar; Cvar; -PHIvar];
J(Dim.Node(3) + 1 : Dim.Node(6), :) = [-Gvar', Cvar', -PHIvar', Hessian + nu_H*eye(Dim.Z)];

KKT_Matrix.J = J;

end

