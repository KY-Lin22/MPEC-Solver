classdef stabilized_SQP_Izamailov2015 < handle
    %Implementation of a globalized and stabilized SQP method which solves 
    %the following degenerate NLP problem:
    %  min f(x),
    %  s.t. h(x) = 0,
    %  g(x)<= 0,
    % Note 1: x \in R^n, f \in R, h \in R^l, g \in R^m
    % Note 2: h and g may be degenerate, i.e., the multipliers associated to
    % stationary points are not unique
    % Ref: A.F.Izmailov et al. "Combining stabilized SQP with the augmented
    % Lagrangian algorithm", Computational Optimization and Applications, 2015£¬
    
    properties
        NLP % struct, symbolic representation and dimmension of the NLP, 
            % with field 'x', 'f', 'h', 'g', 
            %            'fx', 'hx', 'gx', 
            %            'xDim', 'hDim', 'gDim'
            %            'lambda' (Lagrangian multiplier for h), 'mu'(Lagrangian multiplier for g)
            %            'sigma' (dual stabilization parameter), 'eta'(primal direction)
            %            'L', 'Lx', 'Lxx', (Lagrangian function, Jacobian and Hessian)
            %            'AugL', 'AugLx', 'AugLx_eta' (Augmented Lagrangian function, Jacobian and Jacobian vector product)
            %            'rho'(natural residual of KKT system),  'psi'(L2 norm constraint violation),
            %            'qp_H', 'qp_A'(Hessian and Jacobian matrix used in qpoases solver)
        Option % struct, solver option
        FunObj % % struct, CasADi function object 
               % with field 'f', 'h', 'g'
               %            'fx', 'hx', 'gx'
               %            'L', 'Lx', 'Lxx'
               %            'AugL', 'AugLx', 'AugLx_eta'
               %            'rho', 'psi',
               %            'qp_H', 'qp_A'
               %            'qp_solver'
    end
    %% Constructor Method for stabilized_SQP_Izamailov2015
    methods
        function self = stabilized_SQP_Izamailov2015(NLP)
            %stabilized_SQP_Izamailov2015: Construct an instance of this class
            %   Detailed explanation goes here
            
            % import CasADi to workspace
            addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
            import casadi.*        
            
            %% initialize properties: NLP
            % symbolic representation of problem variable, function and Jacobian
            self.NLP = struct('x', NLP.x, 'f', NLP.f, 'h', NLP.h, 'g', NLP.g);
            
            self.NLP.fx = jacobian(self.NLP.f, self.NLP.x);
            self.NLP.hx = jacobian(self.NLP.h, self.NLP.x);    
            self.NLP.gx = jacobian(self.NLP.g, self.NLP.x);
            
            % specify problem dimension 
            self.NLP.xDim = size(self.NLP.x, 1);
            self.NLP.hDim = size(self.NLP.h, 1);
            self.NLP.gDim = size(self.NLP.g, 1);            
            
            % symbolic representation of Lagrangian multiplier, dual stabilization parameter and primal direction eta
            self.NLP.lambda = SX.sym('lambda', self.NLP.hDim, 1);
            self.NLP.mu = SX.sym('mu', self.NLP.gDim, 1);
            self.NLP.sigma = SX.sym('sigma', 1, 1);
            self.NLP.eta = SX.sym('eta', self.NLP.xDim, 1);
            
            % symbolic representation of Lagrangian function, Jacobian and Hessian
            self.NLP.L = self.NLP.f + self.NLP.lambda'*self.NLP.h + self.NLP.mu'*self.NLP.g;                      
            self.NLP.Lx = jacobian(self.NLP.L, self.NLP.x);
            [Lxx, ~] = hessian(self.NLP.L, self.NLP.x);           
            self.NLP.Lxx = Lxx;
            
            % symbolic representation of Augmented Lagrangian function and Jacobian
            self.NLP.AugL = self.NLP.f + (self.NLP.sigma/2)*...
                (norm(self.NLP.lambda + 1/self.NLP.sigma*self.NLP.h)^2 ...
                + norm(max(zeros(self.NLP.gDim, 1), self.NLP.mu + 1/self.NLP.sigma*self.NLP.g))^2);           
            self.NLP.AugLx = jacobian(self.NLP.AugL, self.NLP.x);
            self.NLP.AugLx_eta = jtimes(self.NLP.AugL, self.NLP.x, self.NLP.eta);
            
            % symbolic representation of constraint violation psi and natural residual rho
            self.NLP.psi = norm([self.NLP.h; min(self.NLP.mu ,-self.NLP.g)]);
            self.NLP.rho = norm(self.NLP.Lx) + self.NLP.psi;
            
            % symbolic representation of matrix used in qpoases solver (sparse)
            qp_H = [self.NLP.Lxx,                                          SX(self.NLP.xDim, self.NLP.hDim + self.NLP.gDim);...
                   SX(self.NLP.hDim + self.NLP.gDim, self.NLP.xDim), self.NLP.sigma * SX.eye(self.NLP.hDim + self.NLP.gDim)];           
            self.NLP.qp_H = qp_H;           
            qp_A = [[self.NLP.hx; self.NLP.gx], -self.NLP.sigma * SX.eye(self.NLP.hDim + self.NLP.gDim)];
            self.NLP.qp_A = qp_A;
            
            % display sparsity pattern
            disp('qp_H sparsity pattern: ')
            self.NLP.qp_H.sparsity().spy();           
            disp('qp_A sparsity pattern: ')
            self.NLP.qp_A.sparsity().spy();
            
            %% initialize properties: Option
            self.Option = self.createOption();
            
            %% initialize properties: FunObj
            % function, Jacobian and Hessian for problem and Lagrangian
            self.FunObj.f = Function('f', {self.NLP.x}, {self.NLP.f}, {'x'}, {'f'});
            self.FunObj.h = Function('h', {self.NLP.x}, {self.NLP.h}, {'x'}, {'h'});
            self.FunObj.g = Function('g', {self.NLP.x}, {self.NLP.g}, {'x'}, {'g'});
            
            self.FunObj.fx = Function('fx', {self.NLP.x}, {self.NLP.fx}, {'x'}, {'fx'});
            self.FunObj.hx = Function('hx', {self.NLP.x}, {self.NLP.hx}, {'x'}, {'hx'});
            self.FunObj.gx = Function('gx', {self.NLP.x}, {self.NLP.gx}, {'x'}, {'gx'});
            
            self.FunObj.L = Function('L', {self.NLP.x, self.NLP.lambda, self.NLP.mu}, {self.NLP.L},...
                {'x', 'lambda', 'mu'}, {'L'});            
            self.FunObj.Lx = Function('Lx', {self.NLP.x, self.NLP.lambda, self.NLP.mu}, {self.NLP.Lx},...
                {'x', 'lambda', 'mu'}, {'Lx'});                    
            self.FunObj.Lxx = Function('Lxx', {self.NLP.x, self.NLP.lambda, self.NLP.mu}, {self.NLP.Lxx},...
                {'x', 'lambda', 'mu'}, {'Lxx'});      
            
            self.FunObj.AugL = Function('AugL', {self.NLP.x, self.NLP.lambda, self.NLP.mu, self.NLP.sigma}, {self.NLP.AugL},...
                {'x', 'lambda', 'mu', 'sigma'}, {'AugL'});
            self.FunObj.AugLx = Function('AugLx', {self.NLP.x, self.NLP.lambda, self.NLP.mu, self.NLP.sigma}, {self.NLP.AugLx},...
                {'x', 'lambda', 'mu', 'sigma'}, {'AugLx'});            
            self.FunObj.AugLx_eta = Function('AugLx_eta', {self.NLP.x, self.NLP.lambda, self.NLP.mu, self.NLP.sigma, self.NLP.eta}, {self.NLP.AugLx_eta},...
                {'x', 'lambda', 'mu', 'sigma', 'eta'}, {'AugLx_eta'});
                      
            self.FunObj.psi = Function('psi', {self.NLP.x, self.NLP.mu}, {self.NLP.psi},...
                {'x', 'mu'}, {'psi'});
            self.FunObj.rho = Function('rho', {self.NLP.x, self.NLP.lambda, self.NLP.mu}, {self.NLP.rho},...
                {'x', 'lambda', 'mu'}, {'rho'});
            
            % matrix for qpoases solver (sparse)
            self.FunObj.qp_H = Function('qp_H', {self.NLP.x, self.NLP.lambda, self.NLP.mu, self.NLP.sigma}, {self.NLP.qp_H},...
                {'x', 'lambda', 'mu', 'sigma'}, {'qp_H'});
            self.FunObj.qp_A = Function('qp_A', {self.NLP.x, self.NLP.sigma}, {self.NLP.qp_A},...
                {'x', 'sigma'}, {'qp_A'});
            
            % QP solver
            qp = struct('h', self.NLP.qp_H.sparsity(), 'a', self.NLP.qp_A.sparsity());
            qp_options = struct();     
            qp_options.error_on_fail = false;
            qp_options.printLevel = 'none'; % 'none', 'low', 'medium', 'high' (see qpoases manual sec 5.2)
            qp_options.hessian_type = 'indef';% 'unknown', 'posdef', 'semidef', 'indef', 'zero', 'identity' (see qpoases manual sec 4.4, 4.5)
%             qp_options.enableInertiaCorrection = true;

%             qp_options.enableRegularisation = true; %(this three option corresponds to unbounded QP, see qpoases manual sec 4.5)            
%             qp_options.enableFlippingBounds = true;          
%             qp_options.enableFarBounds = true;
            
            self.FunObj.qp_solver = casadi.conic('qp_solver', 'qpoases', qp, qp_options);            
            
        end
        
    end
    
    
    %% Other Methods for stabilized_SQP_Izamailov2015 
    methods
        function Option = createOption(self)
            Option = struct;
            
            Option.r0 = 1e4; % r0 > 0, init value of natural residual tolerance, used and updated in step 2, 7,
                             % default: 1e4, paragraph 2 in page 424
            Option.epsilon0 = 1e2; % epsilon0 > 0, init value of augmented Lagrangian stationary tolerance, used in step 6, updated in step 2,7, 
                                   % default: 1e2, paragraph 2 in page 424
            Option.sigma0 = 1e-4; % sigma0 > 0, init value of dual stabilization parameter, used in step 1(QP), 3,5,6,7(AugL), updated in step 2,7
                                   % default: 1e-4, paragraph 3 in page 424
            Option.gamma = 1; % gamma > 0, weight parameter in equ(11) checking the descent for AugL, used in the step 3 
                              % default: 1, paragraph 3 in page 424
            
            Option.q = 0.5; % 0 < q < 1, weight parameter in equ(10) updating natural residual tolerance rk, used in step 2, 7
                            % default: 0.5, paragraph 2 in page 424
            Option.xita = 0.5; % 0 < xita < 1, weight parameter in step 7 updating augmented Lagrangian stationary tolerance epsilonk, used in step 7
                              % default: 0.5, paragraph 2 in page 424
            Option.tau = 0.5; % 0 < tau < 1, weight parameter for step size in line search, used in step 5
                             % default: 0.5,  paragraph 3 in page 424
            Option.epsilon = 0.1; % 0 < epsilon < 1, weight parameter for merit function in line search, used in step 5
                                  % default: 0.1,  paragraph 3 in page 424
            Option.kappa = 0.1; % 0 < kappa < 1, weight parameter in equ(15) for updating dual stabilization parameter sigma_k, used in step 7
                                % default: 0.1,  paragraph 3 in page 424
            Option.delta = 0.5; % 0 < delta < 1, weight parameter in equ(15) for updating dual stabilization parameter sigma_k, used in step 7
                                % default: 0.5,  paragraph 3 in page 424
            
            % bounds to safeguard the dual iterate
            Option.bar_lambdaMin = -1e10; % paragraph 3 in page 424
            Option.bar_lambdaMax = 1e10; % paragraph 3 in page 424
            Option.bar_muMax = 1e10; % paragraph 3 in page 424
            
            % modified Hessian
            Option.omegaMax = 1e10;
            
            % termination condition
            Option.maxOuterIterNum = 500; % kMax, positive int(default: 500, end of paragraph 1 in page 424)
            Option.maxInnerIterNum = 50; % jMax, positive int(default: none, so I specify by myself, although 
                                         % proposition 1 shows that j must be finite)
            Option.maxLineSearchIterNum = 8; %iMax, nonnegative int(default: none, so I specify by myself)
            Option.tol_rho = 1e-6; % tolerance for natural residual (default: 1e-6, end of paragraph 1 in page 424)
        end               
        
        function [hat_x_kNext, bar_lambda_kNext, bar_mu_kNext, r_kNext, epsilon_kNext, sigma_kNext, iterateType_kNext, innerLoopFlag] = ...
                innerLoopIteration(self, hat_x_k, bar_lambda_k, bar_mu_k, r_k, epsilon_k, sigma_k, iterateType_k)          
                 
            % x_j: given iterate x_{j}, x_jNext: new iterate x_{j+1}
            % initialize given iterate used in inner loop iteration
            hat_x_j = hat_x_k;            
            has_find_new_outer_iterate = false;% step 2 and 7 will set it as true
            lineSearchFlag = true;
            modifiedHessianFlag = true;
            
            for j = 0 : self.Option.maxInnerIterNum
                %% Checking termination and exitFlag
                % check inner loop termination
                if has_find_new_outer_iterate
                    exitFlag = true;
                    innerLoopFlag = 'Success';
                elseif (j == self.Option.maxInnerIterNum) || (~lineSearchFlag) || (~modifiedHessianFlag)
                    exitFlag = true;
                    innerLoopFlag = 'Fail';
                    disp(' ')
                    if j == self.Option.maxInnerIterNum
                        disp('MSG: inner loop fails because it has reached the maximum inner iteration number')
                    end
                    if ~lineSearchFlag
                        disp('MSG: inner loop fails because line search can not find a not so small step size')
                    end
                    if~modifiedHessianFlag
                        disp('MSG: inner loop fails because omega in modified Hessian exceeds the limit omegaMax')
                    end
                else
                    exitFlag = false;
                end
                % check exitFlag
                if exitFlag
                    % return new outer loop iterate and then break inner loop iteration
                    switch innerLoopFlag
                        case 'Success'
                            % return iterate estimated from step 2 or 7
                            hat_x_kNext = hat_x_kNext_est;
                            bar_lambda_kNext = bar_lambda_kNext_est;
                            bar_mu_kNext = bar_mu_kNext_est;
                            r_kNext = r_kNext_est;
                            epsilon_kNext = epsilon_kNext_est;
                            sigma_kNext = sigma_kNext_est;
                            iterateType_kNext = iterateType_kNext_est;
                        case 'Fail'
                            % return iterate from input
                            hat_x_kNext = hat_x_k;
                            bar_lambda_kNext = bar_lambda_k;
                            bar_mu_kNext = bar_mu_k;
                            r_kNext = r_k;
                            epsilon_kNext = epsilon_k;
                            sigma_kNext = sigma_k;
                            iterateType_kNext = iterateType_k;
                    end
                    break
                end              
                
                %% function and Jacobian evaluation of given inner loop iterate
                % problem function and Jacobian
                fx_j = self.FunObj.fx(hat_x_j);
                h_j  = self.FunObj.h(hat_x_j);
                g_j  = self.FunObj.g(hat_x_j);
                % matrix used in qpoases
                qp_H_j = self.FunObj.qp_H(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k);
                qp_A_j = self.FunObj.qp_A(hat_x_j, sigma_k);
                qp_g_j = [fx_j, casadi.DM.zeros(1, self.NLP.hDim + self.NLP.gDim)];
                qp_lba_j = [-h_j - sigma_k * bar_lambda_k;...
                    -inf*casadi.DM.ones(self.NLP.gDim, 1)];
                qp_uba_j = [-h_j - sigma_k * bar_lambda_k;...
                    -g_j - sigma_k * bar_mu_k];
                HessianType = 'LagrangianHessian';
                
                %% Algorithm 1 (step 1 - 4 while loop)
                % break when it has found a new iterate (step 2) or a direction which is descent for AugL (step 3)
                flag_find_new_outer_iterate_or_descent_for_AugL = false;
                omega = 10;
                while ~flag_find_new_outer_iterate_or_descent_for_AugL
                    % step 1: solve QP subproblem equ (7) for primal and dual direction
                    qp_solution = self.FunObj.qp_solver('h', qp_H_j, 'g', qp_g_j, 'a', qp_A_j, ...
                        'lba', qp_lba_j, 'uba', qp_uba_j,...
                        'x0', zeros(self.NLP.xDim + self.NLP.hDim + self.NLP.gDim, 1));
                    dx = qp_solution.x(1 : self.NLP.xDim, 1);
                    dlambda = qp_solution.x(self.NLP.xDim + 1 : self.NLP.xDim + self.NLP.hDim, 1) - bar_lambda_k;
                    dmu = qp_solution.x(self.NLP.xDim + self.NLP.hDim + 1 : end, 1) - bar_mu_k;
                    
                    qp_status = self.FunObj.qp_solver.stats.return_status;
                    disp(['inner iter j = ', num2str(j), ' --> qpoases return status: ', qp_status])
                    if strcmp(qp_status, 'Successful return.')
                        % step 2: check whether we can obtain a new iterate with full step along the direction
                        check_lambda =  (full(min(bar_lambda_k + dlambda)) >= self.Option.bar_lambdaMin) && ...
                            (full(max(bar_lambda_k + dlambda)) <= self.Option.bar_lambdaMax) ;
                        check_mu = (full(min(bar_mu_k + dmu)) >= 0) && (full(max(bar_mu_k + dmu)) <= self.Option.bar_muMax) ;
                        rho_est = self.FunObj.rho(hat_x_j + dx, bar_lambda_k + dlambda, bar_mu_k + dmu);
                        equ8_is_satisfied = (check_lambda) && (check_mu) && (full(rho_est) <= full(r_k));
                        if (strcmp(HessianType, 'LagrangianHessian')) && (equ8_is_satisfied)
                            % set new outer iterate (sSQP iteration) and then break this while loop
                            has_find_new_outer_iterate = true;
                            hat_x_kNext_est = hat_x_j + dx;
                            bar_lambda_kNext_est = bar_lambda_k + dlambda;
                            bar_mu_kNext_est = bar_mu_k + dmu;
                            r_kNext_est = rho_est;
                            epsilon_kNext_est = epsilon_k;
                            sigma_kNext_est = self.Option.q * r_k;    
                            iterateType_kNext_est = 'sSQP_iterate';
                            break
                        end
                        % step 3: check whether we can obtain a primal direction which is descent for augmented Lagrangian
                        AugLx_eta = self.FunObj.AugLx_eta(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k, dx);
                        norm_dx = norm(dx);
                        equ11_is_satisfied = (full(AugLx_eta) <= full(-self.Option.gamma*norm_dx^2));
                        if equ11_is_satisfied
                            % break this while loop and perform step 5
                            break
                        end
                        % step 4: modify Hessian    
                        if omega >= self.Option.omegaMax
                            % break this while loop and terminate the overall rountie
                            modifiedHessianFlag = false;
                            break
                        end
                        modify_Hessian = [omega*casadi.DM.eye(self.NLP.xDim), casadi.DM.zeros(self.NLP.xDim, self.NLP.hDim + self.NLP.gDim);...
                            casadi.DM.zeros(self.NLP.hDim + self.NLP.gDim, self.NLP.xDim), casadi.DM.zeros(self.NLP.hDim + self.NLP.gDim, self.NLP.hDim + self.NLP.gDim)];
                        qp_H_j = qp_H_j + modify_Hessian;
                        HessianType = 'ModifiedHessian';
                        omega = 10 * omega;
                    else
                        % step 4: modify Hessian
                        if omega >= self.Option.omegaMax
                            % break this while loop and terminate the overall rountie
                            modifiedHessianFlag = false;
                            break
                        end
                        modify_Hessian = [omega*casadi.DM.eye(self.NLP.xDim), casadi.DM.zeros(self.NLP.xDim, self.NLP.hDim + self.NLP.gDim);...
                            casadi.DM.zeros(self.NLP.hDim + self.NLP.gDim, self.NLP.xDim), casadi.DM.zeros(self.NLP.hDim + self.NLP.gDim, self.NLP.hDim + self.NLP.gDim)];
                        qp_H_j = qp_H_j + modify_Hessian;
                        HessianType = 'ModifiedHessian';
                        omega = 10 * omega;                        
                    end
                end
                                               
                 %% Algorithm 1 (step 5 - step 7 if loop)
                 % perform AugL iteration after step 3, i.e., the case that step 2 does not satisfy
                 if (~has_find_new_outer_iterate) && modifiedHessianFlag
                     % step 5: line search
                     for i = 0 : self.Option.maxLineSearchIterNum
                         % check failure flag
                         if i == self.Option.maxLineSearchIterNum
                             % break line search 'for loop'
                             lineSearchFlag = false;
                             break
                         end
                         % Amijo rules
                         AugL_j = self.FunObj.AugL(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k);
                         hat_x_jNext = hat_x_j + (self.Option.tau)^i*dx;
                         AugL_jNext = self.FunObj.AugL(hat_x_jNext, bar_lambda_k, bar_mu_k, sigma_k);
                         if (full(AugL_jNext) <= full(AugL_j + self.Option.epsilon * (self.Option.tau)^i * AugLx_eta))
                             % terminate line search 'for loop'
                             break
                         end
                     end
                     
                     if lineSearchFlag
                         % step 6 check stationary of AugL
                         AugLx_jNext = self.FunObj.AugLx(hat_x_jNext, bar_lambda_k, bar_mu_k, sigma_k);
                         if full(norm(AugLx_jNext)) <= epsilon_k
                             % step 7: set new outer iterate (AugL iteration)                             
                             has_find_new_outer_iterate = true;
                             epsilon_kNext_est = self.Option.xita * epsilon_k;
                             hat_x_kNext_est = hat_x_jNext;                            
                             % dual variable
                             h_jNext = self.FunObj.h(hat_x_jNext);
                             bar_lambda_kNext_est = bar_lambda_k + 1/sigma_k*h_jNext;
                             if (full(min(bar_lambda_kNext_est)) >= self.Option.bar_lambdaMin) && (full(max(bar_lambda_kNext_est)) <= self.Option.bar_lambdaMax)
                                 % bar_lambda_kNext_est satisfies dual bounds
                             else
                                 bar_lambda_kNext_est = min(max(self.Option.bar_lambdaMin*DM.ones(self.NLP.hDim, 1), bar_lambda_kNext_est),...
                                     self.Option.bar_lambdaMax *DM.ones(self.NLP.hDim, 1));
                             end
                             g_jNext = self.FunObj.g(hat_x_jNext);
                             bar_mu_kNext_est = max(zeros(self.NLP.gDim, 1), bar_mu_k + 1/sigma_k*g_jNext);
                             if (full(min(bar_mu_kNext_est)) >= 0) && (full(max(bar_mu_kNext_est)) <= self.Option.bar_muMax)
                                 % bar_mu_kNext_est satisfies dual bounds
                             else
                                 bar_mu_kNext_est = min(max(DM.zeros(self.NLP.gDim, 1), bar_mu_kNext_est),...
                                     self.Option.bar_muMax *DM.ones(self.NLP.gDim, 1));
                             end
                             % natural residual tol and dual stabilization variable
                             rho_kNext_est = self.FunObj.rho(hat_x_kNext_est, bar_lambda_kNext_est, bar_mu_kNext_est);
                             if full(rho_kNext_est) <= full(r_k)
                                 r_kNext_est = self.Option.q * r_k;
                                 sigma_kNext_est = rho_kNext_est;
                             else
                                 r_kNext_est = r_k;
                                 psi_kNext_est = self.FunObj.psi(hat_x_kNext_est, bar_mu_kNext_est);
                                 psi_k = self.FunObj.psi(hat_x_k, bar_mu_k);
                                 if full(psi_kNext_est) <= full(self.Option.delta*psi_k)
                                     sigma_kNext_est = sigma_k;
                                 else
                                     sigma_kNext_est = self.Option.kappa*sigma_k;
                                 end
                             end
                             iterateType_kNext_est = 'AugL_iterate';
                         else
                             % prepare primal variable for next inner loop iteration
                             hat_x_j = hat_x_jNext;
                         end
                         
                     end
                     
                 end
                           
            end
            
        end
        
        function [x_Opt, Info] = solveNLP(self, x_Init)
            %%                  
            timeStart = tic;
            % very first primal and dual variable
            x_0 = x_Init;
            lambda_0 = zeros(self.NLP.hDim, 1); % here I initilize the very first dual variables with value which satisfies their estimate    
            mu_0 = zeros(self.NLP.gDim, 1);     % bounds (paragraph 3 in page 424) so that their first estimate can choose these value. 
            
            % initialize estimate for the first dual variable (here 0 is the value of outer loop counter k)
            bar_lambda_0 = lambda_0; % \in [Option.bar_lambdaMin, Option.bar_lambdaMax]
            bar_mu_0 = mu_0; % \in [0, Option.bar_muMax]
            
            % initialize estimate for the first primal variable (here 0 is the value of outer loop counter k and inner loop counter j)
            hat_x_0 = x_0;
            
            %% SQP iteration rountie
            % k: outer iteration counter, increase after sSQP(step 2 in ref.) and Aug-L iteration(step 6 and 7 in ref.)
            % j: inner iteration counter, increase after current Aug-L subproblem does not have acceptable stationary property(step 6 in ref.)            
            % in k iteration, x_k is the given iterate x_{k}, x_kNext is the new iterate x_{k+1}
            
            % initialize given iterate used in outer loop iteration
            hat_x_k = hat_x_0;
            bar_lambda_k = bar_lambda_0;
            bar_mu_k = bar_mu_0;
            r_k = self.Option.r0;
            epsilon_k = self.Option.epsilon0;
            sigma_k = self.Option.sigma0;
            iterateType_k = 'Init';
            innerLoopFlag = 'Success';
            disp('********************************************************************')
            % outer loop iteration
            for k = 0 : self.Option.maxOuterIterNum
                % STEP 1: compute natural residual rho of given iterate (k)
                rho_k = self.FunObj.rho(hat_x_k, bar_lambda_k, bar_mu_k);                
                % STEP 2: check outer loop termination (end of paragraph 1 in page 424)               
                disp(['Iter: ', num2str(k), '; ',...
                    'iterateType: ', iterateType_k, '; ',...
                    'rho: ', num2str(full(rho_k))])
                disp('--------------------------------------------------------------------')
                if full(rho_k) < self.Option.tol_rho
                    % solver finds the optimal solution
                    exitFlag = true;
                    terminalStatus = 'Success';
                elseif (k == self.Option.maxOuterIterNum) || (strcmp(innerLoopFlag, 'Fail'))
                    % solver fails to find the optimal solution
                    exitFlag = true;
                    terminalStatus = 'Fail'; 
                else
                    exitFlag = false;
                end
                % STEP 3: check exitFlag
                if exitFlag
                    timeElapsed = toc(timeStart);
                    % return given iterate as solution
                    x_Opt = hat_x_k;
                    % organize output information   
                    Info = struct('iterNum', k, 'terminalStatus', terminalStatus, 'time', timeElapsed, ...
                        'lambda', bar_lambda_k, 'mu', bar_mu_k,...
                        'rho', rho_k, 'iterateType', iterateType_k);
                    % print information
                    disp('********************************************************************')
                    disp('Done!')
                    disp(['iterNum: ', num2str(Info.iterNum)])
                    disp(['terminalStatus: ', Info.terminalStatus])                    
                    disp(['elapsedTime: ', num2str(Info.time), ' seconds'])
                    disp(['rho: ', num2str(full(Info.rho))])
                    disp(['solution type: ', Info.iterateType])
                    break
                end
                % STEP 4: inner loop iteration (Algorithm 1)
                [hat_x_kNext, bar_lambda_kNext, bar_mu_kNext, r_kNext, epsilon_kNext, sigma_kNext, iterateType_kNext, innerLoopFlag] = ...
                    self.innerLoopIteration(hat_x_k, bar_lambda_k, bar_mu_k, r_k, epsilon_k, sigma_k, iterateType_k);
                % STEP 5: prepare for next outer iteration
                hat_x_k = hat_x_kNext;
                bar_lambda_k = bar_lambda_kNext;
                bar_mu_k = bar_mu_kNext;
                r_k = r_kNext;
                epsilon_k = epsilon_kNext;
                sigma_k = sigma_kNext;
                iterateType_k = iterateType_kNext;
            end
            
        end
               
    end
end

