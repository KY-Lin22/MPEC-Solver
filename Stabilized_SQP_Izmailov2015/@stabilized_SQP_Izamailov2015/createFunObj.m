function FunObj = createFunObj(self)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
import casadi.* 
NLP = self.NLP;
% function, Jacobian and Hessian for problem and Lagrangian
FunObj.f = Function('f', {NLP.x}, {NLP.f}, {'x'}, {'f'});
FunObj.h = Function('h', {NLP.x}, {NLP.h}, {'x'}, {'h'});
FunObj.g = Function('g', {NLP.x}, {NLP.g}, {'x'}, {'g'});

FunObj.fx = Function('fx', {NLP.x}, {NLP.fx}, {'x'}, {'fx'});
FunObj.hx = Function('hx', {NLP.x}, {NLP.hx}, {'x'}, {'hx'});
FunObj.gx = Function('gx', {NLP.x}, {NLP.gx}, {'x'}, {'gx'});

FunObj.L = Function('L', {NLP.x, NLP.lambda, NLP.mu}, {NLP.L},...
    {'x', 'lambda', 'mu'}, {'L'});
FunObj.Lx = Function('Lx', {NLP.x, NLP.lambda, NLP.mu}, {NLP.Lx},...
    {'x', 'lambda', 'mu'}, {'Lx'});
FunObj.Lxx = Function('Lxx', {NLP.x, NLP.lambda, NLP.mu}, {NLP.Lxx},...
    {'x', 'lambda', 'mu'}, {'Lxx'});

FunObj.AugL = Function('AugL', {NLP.x, NLP.lambda, NLP.mu, NLP.sigma}, {NLP.AugL},...
    {'x', 'lambda', 'mu', 'sigma'}, {'AugL'});
FunObj.AugLx = Function('AugLx', {NLP.x, NLP.lambda, NLP.mu, NLP.sigma}, {NLP.AugLx},...
    {'x', 'lambda', 'mu', 'sigma'}, {'AugLx'});
FunObj.AugLx_eta = Function('AugLx_eta', {NLP.x, NLP.lambda, NLP.mu, NLP.sigma, NLP.eta}, {NLP.AugLx_eta},...
    {'x', 'lambda', 'mu', 'sigma', 'eta'}, {'AugLx_eta'});

FunObj.psi = Function('psi', {NLP.x, NLP.mu}, {NLP.psi},...
    {'x', 'mu'}, {'psi'});
FunObj.rho = Function('rho', {NLP.x, NLP.lambda, NLP.mu}, {NLP.rho},...
    {'x', 'lambda', 'mu'}, {'rho'});

% matrix for qpoases solver (sparse)
FunObj.qp_H = Function('qp_H', {NLP.x, NLP.lambda, NLP.mu, NLP.sigma}, {NLP.qp_H},...
    {'x', 'lambda', 'mu', 'sigma'}, {'qp_H'});
FunObj.qp_A = Function('qp_A', {NLP.x, NLP.sigma}, {NLP.qp_A},...
    {'x', 'sigma'}, {'qp_A'});

% QP solver
qp = struct('h', NLP.qp_H.sparsity(), 'a', NLP.qp_A.sparsity());
qp_options = struct();
qp_options.error_on_fail = false;
qp_options.printLevel = 'none'; % 'none', 'low', 'medium', 'high' (see qpoases manual sec 5.2)
qp_options.hessian_type = 'indef';% 'unknown', 'posdef', 'semidef', 'indef', 'zero', 'identity' (see qpoases manual sec 4.4, 4.5)
% qp_options.enableInertiaCorrection = true;

% qp_options.enableRegularisation = true; %(this three option corresponds to unbounded QP, see qpoases manual sec 4.5)
% qp_options.enableFlippingBounds = true;
% qp_options.enableFarBounds = true;

FunObj.qp_solver = conic('qp_solver', 'qpoases', qp, qp_options);
end

