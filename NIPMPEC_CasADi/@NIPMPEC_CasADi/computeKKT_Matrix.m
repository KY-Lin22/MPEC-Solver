function KKT_Matrix = computeKKT_Matrix(self, Jac, Hessian)
%UNTITLED30 Summary of this function goes here
%   Detailed explanation goes here

Dim = self.MPEC.Dim;
% regular parameter
nu_G = self.Option.RegularParam.nu_G;
nu_J = self.Option.RegularParam.nu_J;
nu_H = self.Option.RegularParam.nu_H;

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

if strcmp(self.Option.linearSystemSolver, 'mldivide_sparse')
    %% sparse
    % diagVec
    i_diagVec = 1 : Dim.Node(3);
    j_diagVec = i_diagVec;
    s_diagVec = diagVec;
    % [-Gx, -Gp, -Gw]
    [i_Gx, j_Gx, s_Gx] = find(-Jac.Gx);
    j_Gx = j_Gx + Dim.Node(3);
    [i_Gp, j_Gp, s_Gp] = find(-Jac.Gp);
    j_Gp = j_Gp + Dim.Node(4);
    [i_Gw, j_Gw, s_Gw] = find(-Jac.Gw);
    j_Gw = j_Gw + Dim.Node(5);    
    % [Cx, Cp, Cw]
    [i_Cx, j_Cx, s_Cx] = find(Jac.Cx);
    i_Cx = i_Cx + Dim.Node(1);
    j_Cx = j_Cx + Dim.Node(3);
    [i_Cp, j_Cp, s_Cp] = find(Jac.Cp);
    i_Cp = i_Cp + Dim.Node(1);
    j_Cp = j_Cp + Dim.Node(4);  
    [i_Cw, j_Cw, s_Cw] = find(Jac.Cw);
    i_Cw = i_Cw + Dim.Node(1);
    j_Cw = j_Cw + Dim.Node(5); 
    % [-PHIx, -PHIp, -PHIw]
    [i_PHIx, j_PHIx, s_PHIx] = find(-Jac.PHIx);
    i_PHIx = i_PHIx + Dim.Node(2);
    j_PHIx = j_PHIx + Dim.Node(3);
    [i_PHIp, j_PHIp, s_PHIp] = find(-Jac.PHIp);
    i_PHIp = i_PHIp + Dim.Node(2);
    j_PHIp = j_PHIp + Dim.Node(4);
    [i_PHIw, j_PHIw, s_PHIw] = find(-Jac.PHIw);
    i_PHIw = i_PHIw + Dim.Node(2);
    j_PHIw = j_PHIw + Dim.Node(5);   
       
    % [-Gx'; -Gp'; -Gw']
    i_GxT = j_Gx;
    j_GxT = i_Gx;
    s_GxT = s_Gx;
    i_GpT = j_Gp;   
    j_GpT = i_Gp;
    s_GpT = s_Gp; 
    i_GwT = j_Gw;  
    j_GwT = i_Gw;
    s_GwT = s_Gw;
    % [Cx'; Cp'; Cw']
    i_CxT = j_Cx;
    j_CxT = i_Cx;
    s_CxT = s_Cx;
    i_CpT = j_Cp;
    j_CpT = i_Cp;
    s_CpT = s_Cp;
    i_CwT = j_Cw;
    j_CwT = i_Cw; 
    s_CwT = s_Cw;
    % [-PHIx'; -PHIp'; -PHIw']
    i_PHIxT = j_PHIx;
    j_PHIxT = i_PHIx;
    s_PHIxT = s_PHIx;
    i_PHIpT = j_PHIp;
    j_PHIpT = i_PHIp;   
    s_PHIpT = s_PHIp;
    i_PHIwT = j_PHIw;
    j_PHIwT = i_PHIw; 
    s_PHIwT = s_PHIw;
    
    % Hessian
    Hessian_diag = diag(Hessian);
    i_H_diag = Dim.Node(3) + 1 : Dim.Node(6);
    j_H_diag = i_H_diag;
    s_H_diag = Hessian_diag + nu_H*ones(Dim.Z, 1);
    
    Hessian_tril = tril(Hessian, -1);
    [i_H_tril, j_H_tril, s_H_tril] = find(Hessian_tril);
    i_H_tril = i_H_tril + Dim.Node(3);
    j_H_tril = j_H_tril + Dim.Node(3);
    
    i_H_tril_T = j_H_tril;
    j_H_tril_T = i_H_tril;
    s_H_tril_T = s_H_tril;    

    % J
    i = [i_diagVec';...
        i_Gx;  i_Gp;  i_Gw;  i_Cx;  i_Cp;  i_Cw;  i_PHIx;  i_PHIp;  i_PHIw;...
        i_GxT; i_GpT; i_GwT; i_CxT; i_CpT; i_CwT; i_PHIxT; i_PHIpT; i_PHIwT;...
        i_H_diag'; i_H_tril; i_H_tril_T];
    j = [j_diagVec';...
        j_Gx;  j_Gp;  j_Gw;  j_Cx;  j_Cp;  j_Cw;  j_PHIx;  j_PHIp;  j_PHIw;...
        j_GxT; j_GpT; j_GwT; j_CxT; j_CpT; j_CwT; j_PHIxT; j_PHIpT; j_PHIwT;...
        j_H_diag'; j_H_tril; j_H_tril_T];
    s = [s_diagVec;...
        s_Gx;  s_Gp;  s_Gw;  s_Cx;  s_Cp;  s_Cw;  s_PHIx;  s_PHIp;  s_PHIw;...
        s_GxT; s_GpT; s_GwT; s_CxT; s_CpT; s_CwT; s_PHIxT; s_PHIpT; s_PHIwT;...
        s_H_diag; s_H_tril; s_H_tril_T];
    J = sparse(i, j, s, Dim.Y, Dim.Y);

else
    %% dense
    % Constraint Jacobian
    ConstraintJacobian =...
        [-Jac.Gx,   -Jac.Gp,   -Jac.Gw;...
          Jac.Cx,    Jac.Cp,    Jac.Cw;...
         -Jac.PHIx, -Jac.PHIp, -Jac.PHIw];
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

