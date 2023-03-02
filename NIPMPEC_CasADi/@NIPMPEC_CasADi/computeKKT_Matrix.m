function KKT_Matrix = computeKKT_Matrix(self, Jac, Hessian)
%UNTITLED30 Summary of this function goes here
%   Detailed explanation goes here

Dim = self.MPEC.Dim;
Option = self.Option;
% regular parameter
nu_G = Option.RegularParam.nu_G;
nu_J = Option.RegularParam.nu_J;
nu_H = Option.RegularParam.nu_H;

% diag matrix vector: D = diag(d), nu_J,  E = diag(e);
d = -(Jac.PSIgSigma_diagVec - nu_G * ones(Dim.sigma, 1))./(Jac.PSIgG_diagVec - nu_G * ones(Dim.sigma, 1));
nu_J_vec = -nu_J*ones(Dim.eta, 1);
e = -(Jac.PSIphiGamma_diagVec - nu_G * ones(Dim.gamma, 1))./(Jac.PSIphiPHI_diagVec - nu_G * ones(Dim.gamma, 1));
diagVec = [d; nu_J_vec; e];

%% assemble KKT matrix
% J = [D,                           zeros(Dim.sigma, Dim.eta),  zeros(Dim.sigma, Dim.gamma),  -G_grad;...
%     zeros(Dim.eta, Dim.sigma),    -nu_J * eye(Dim.eta),       zeros(Dim.eta, Dim.gamma),    C_grad;...
%     zeros(Dim.gamma, Dim.sigma),  zeros(Dim.gamma, Dim.eta),  E,                            -PHI_grad;...
%     -G_grad',                     C_grad',                    -PHI_grad',                   Hessian + nu_H*eye(Dim.Z)]; 

if strcmp(Option.linearSystemSolver, 'mldivide_sparse')
    %% sparse
    % diagVec
    i_diagVec = 1 : Dim.Node(3);
    j_diagVec = i_diagVec;
    s_diagVec = diagVec;

    % -G_grad = [-Gx, -Gp, ~]
    [i_Gx, j_Gx, s_Gx] = find(-Jac.Gx);
    j_Gx = j_Gx + Dim.Node(3);
    [i_Gp, j_Gp, s_Gp] = find(-Jac.Gp);
    j_Gp = j_Gp + Dim.Node(4);
   
    % C_grad = [Cx, Cp, Cw]
    [i_Cx, j_Cx, s_Cx] = find(Jac.Cx);
    i_Cx = i_Cx + Dim.Node(1);
    j_Cx = j_Cx + Dim.Node(3);
    [i_Cp, j_Cp, s_Cp] = find(Jac.Cp);
    i_Cp = i_Cp + Dim.Node(1);
    j_Cp = j_Cp + Dim.Node(4);  
    [i_Cw, j_Cw, s_Cw] = find(Jac.Cw);
    i_Cw = i_Cw + Dim.Node(1);
    j_Cw = j_Cw + Dim.Node(5); 
    % PHI_grad = [~, -PHIp, -PHIw]
    [i_PHIp, j_PHIp, s_PHIp] = find(-Jac.PHIp);
    i_PHIp = i_PHIp + Dim.Node(2);
    j_PHIp = j_PHIp + Dim.Node(4);
    [i_PHIw, j_PHIw, s_PHIw] = find(-Jac.PHIw);
    i_PHIw = i_PHIw + Dim.Node(2);
    j_PHIw = j_PHIw + Dim.Node(5);   
       
    % -G_grad' = [-Gx'; -Gp'; ~]
    i_GxT = j_Gx;
    j_GxT = i_Gx;
    s_GxT = s_Gx;
    i_GpT = j_Gp;   
    j_GpT = i_Gp;
    s_GpT = s_Gp; 
    % C_grad' = [Cx'; Cp'; Cw']
    i_CxT = j_Cx;
    j_CxT = i_Cx;
    s_CxT = s_Cx;
    i_CpT = j_Cp;
    j_CpT = i_Cp;
    s_CpT = s_Cp;
    i_CwT = j_Cw;
    j_CwT = i_Cw; 
    s_CwT = s_Cw;
    % -PHI_grad' = [~; -PHIp'; -PHIw']
    i_PHIpT = j_PHIp;
    j_PHIpT = i_PHIp;   
    s_PHIpT = s_PHIp;
    i_PHIwT = j_PHIw;
    j_PHIwT = i_PHIw; 
    s_PHIwT = s_PHIw;

    % Hessian
    [i_H, j_H, s_H] = find(Hessian + nu_H*speye(Dim.Z));
    i_H = i_H + Dim.Node(3);
    j_H = j_H + Dim.Node(3);

    % J
    i = [i_diagVec';...
        i_Gx;  i_Gp;  i_Cx;  i_Cp;  i_Cw;  i_PHIp;  i_PHIw;...
        i_GxT; i_GpT; i_CxT; i_CpT; i_CwT; i_PHIpT; i_PHIwT;...
        i_H];
    j = [j_diagVec';...
        j_Gx;  j_Gp;  j_Cx;  j_Cp;  j_Cw;  j_PHIp;  j_PHIw;...
        j_GxT; j_GpT; j_CxT; j_CpT; j_CwT; j_PHIpT; j_PHIwT;...
        j_H];
    s = [s_diagVec;...
        s_Gx;  s_Gp;  s_Cx;  s_Cp;  s_Cw;  s_PHIp;  s_PHIw;...
        s_GxT; s_GpT; s_CxT; s_CpT; s_CwT; s_PHIpT; s_PHIwT;...
        s_H];
    J = sparse(i, j, s, Dim.Y, Dim.Y, length(s));

else
    %% dense
    % Constraint Jacobian
    ConstraintJacobian =...
        [-Jac.Gx,                 -Jac.Gp,    zeros(Dim.sigma, Dim.w);...
          Jac.Cx,                  Jac.Cp,    Jac.Cw;...
         zeros(Dim.gamma, Dim.x), -Jac.PHIp, -Jac.PHIw];
    % J
    J = zeros(Dim.Y, Dim.Y);
    for i = 1 : Dim.Node(3)
        J(i, i) = diagVec(i);
    end
    J(1 : Dim.Node(3), Dim.Node(3) + 1 : Dim.Node(6)) = ConstraintJacobian;
    J(Dim.Node(3) + 1 : Dim.Node(6), :) = [ConstraintJacobian', Hessian + nu_H*eye(Dim.Z)];
end

KKT_Matrix.J = J;

end

