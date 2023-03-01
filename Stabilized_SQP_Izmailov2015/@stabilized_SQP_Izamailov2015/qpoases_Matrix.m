function [qp_H_j, qp_A_j, qp_g_j, qp_lba_j, qp_uba_j] = qpoases_Matrix(self, hat_x_j, bar_lambda_k, bar_mu_k, sigma_k)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
NLP = self.NLP;
FunObj = self.FunObj;

% problem function and Jacobian
fx_j = FunObj.fx(hat_x_j);
h_j  = FunObj.h(hat_x_j);
g_j  = FunObj.g(hat_x_j);

fx_j = full(fx_j);
h_j = full(h_j);
g_j = full(g_j);

% matrix used in qpoases
qp_H_j = FunObj.qp_H(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k);
qp_A_j = FunObj.qp_A(hat_x_j, sigma_k);
qp_H_j = full(qp_H_j);
qp_A_j = full(qp_A_j);

qp_g_j = [fx_j, zeros(1, NLP.hDim + NLP.gDim)];
qp_lba_j = [-h_j - sigma_k * bar_lambda_k;...
    -inf*ones(NLP.gDim, 1)];
qp_uba_j = [-h_j - sigma_k * bar_lambda_k;...
    -g_j - sigma_k * bar_mu_k];

end

