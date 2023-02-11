function [s_k, z_k, Fun_k] = computePerturedParam(self, Var_k, Fun_k, s, z)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
Option = self.Option;
sEnd = Option.sEnd;
zEnd = Option.zEnd;

kappa_s_times = 0.2;
kappa_s_exp = 1.5;
kappa_z_times = 0.2;
kappa_z_exp = 1.5;

% set scaling parameter and threshold 
kappa_F = 10; 
threshold_sz = max([s, z]);

% constraint violation (L Inf norm)
totalCstrVio_LInfNorm = norm([Fun_k.PSIg; Fun_k.C; Fun_k.PSIphi], Inf);

%% update s and z
if totalCstrVio_LInfNorm < kappa_F * threshold_sz
    % close to the optimal solution, need to set them a smaller value
    s_k_trail = min([kappa_s_times .* s, s.^kappa_s_exp]);
    s_k = max([s_k_trail, sEnd]);
    z_k_trail = min([kappa_z_times .* z, z.^kappa_z_exp]);
    z_k = max([z_k_trail, zEnd]);
else
    % far way from the optimal solution, set them the previous value
    s_k = s;
    z_k = z;
end

%% update Fun_k
if (s_k == s) && (z_k == z)
    % both s and z do not update, hence no function need to be updated
    
elseif (s_k ~= s) && (z_k == z)
    % only s is update, hence update function about PHI
    % PHI
    PHI_k = self.FunObj.PHI(Var_k.x, Var_k.p, Var_k.w, s_k);
    Fun_k.PHI = full(PHI_k);
    % FB for PHI
    PSIphi_k = self.FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphi = full(PSIphi_k);      
    [PSIphiGamma_diagVec_k, PSIphiPHI_diagVec_k] = self.FunObj.FB_PHI_grad(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphiGamma = diag(full(PSIphiGamma_diagVec_k));
    Fun_k.PSIphiPHI = diag(full(PSIphiPHI_diagVec_k));    
       
elseif (s_k == s) && (z_k ~= z)
    % only z is update, hence update function about FB
    % FB for G
    PSIg_k = self.FunObj.FB_G(Var_k.sigma, Fun_k.G, z_k);
    Fun_k.PSIg = full(PSIg_k);    
    [PSIgSigma_diagVec_k, PSIgG_k_diagVec] = self.FunObj.FB_G_grad(Var_k.sigma, Fun_k.G, z_k);
    Fun_k.PSIgSigma = diag(full(PSIgSigma_diagVec_k));
    Fun_k.PSIgG = diag(full(PSIgG_k_diagVec));
    % FB for PHI
    PSIphi_k = self.FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphi = full(PSIphi_k);    
    [PSIphiGamma_diagVec_k, PSIphiPHI_diagVec_k] = self.FunObj.FB_PHI_grad(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphiGamma = diag(full(PSIphiGamma_diagVec_k));
    Fun_k.PSIphiPHI = diag(full(PSIphiPHI_diagVec_k));
    
else
    % both s and z update, hence update function about PHI and FB
    % PHI
    PHI_k = self.FunObj.PHI(Var_k.x, Var_k.p, Var_k.w, s_k);
    Fun_k.PHI = full(PHI_k);
    % FB for G
    PSIg_k = self.FunObj.FB_G(Var_k.sigma, Fun_k.G, z_k);
    Fun_k.PSIg = full(PSIg_k);    
    [PSIgSigma_diagVec_k, PSIgG_k_diagVec] = self.FunObj.FB_G_grad(Var_k.sigma, Fun_k.G, z_k);
    Fun_k.PSIgSigma = diag(full(PSIgSigma_diagVec_k));
    Fun_k.PSIgG = diag(full(PSIgG_k_diagVec));
    % FB for PHI
    PSIphi_k = self.FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphi = full(PSIphi_k);    
    [PSIphiGamma_diagVec_k, PSIphiPHI_diagVec_k] = self.FunObj.FB_PHI_grad(Var_k.gamma, Fun_k.PHI, z_k);
    Fun_k.PSIphiGamma = diag(full(PSIphiGamma_diagVec_k));
    Fun_k.PSIphiPHI = diag(full(PSIphiPHI_diagVec_k));    
        
end

end

