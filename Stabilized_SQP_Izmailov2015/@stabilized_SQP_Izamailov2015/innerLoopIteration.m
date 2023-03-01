function [hat_x_kNext, bar_lambda_kNext, bar_mu_kNext, r_kNext, epsilon_kNext, sigma_kNext, iterateType_kNext, innerLoopFlag] = ...
    innerLoopIteration(self, hat_x_k, bar_lambda_k, bar_mu_k, r_k, epsilon_k, sigma_k, iterateType_k)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
NLP = self.NLP;
Option = self.Option;
FunObj = self.FunObj;

% x_j: given iterate x_{j}, x_jNext: new iterate x_{j+1}
% initialize given iterate used in inner loop iteration
hat_x_j = hat_x_k;
has_find_new_outer_iterate = false;% step 2 and 7 will set it as true
lineSearchFlag = true;
modifiedHessianFlag = true;

for j = 0 : Option.maxInnerIterNum
    %% Checking termination
    % check inner loop termination
    if has_find_new_outer_iterate
        exitFlag = true;
        innerLoopFlag = 'Success';
    elseif (j == Option.maxInnerIterNum) || (~lineSearchFlag) || (~modifiedHessianFlag)
        exitFlag = true;
        innerLoopFlag = 'Fail';
        disp(' ')
        if j == Option.maxInnerIterNum
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
    %% check exitFlag
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
    
    %% matrix evaluation of given inner loop iterate for qpoases
    [qp_H_j, qp_A_j, qp_g_j, qp_lba_j, qp_uba_j] = self.qpoases_Matrix(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k);

    HessianType = 'LagrangianHessian';
    
    %% Algorithm 1 (step 1 - 4 while loop)
    % break when it has found a new iterate (step 2) or a direction which is descent for AugL (step 3)
    flag_find_new_outer_iterate_or_descent_for_AugL = false;
    omega = 10;
    while ~flag_find_new_outer_iterate_or_descent_for_AugL
        % step 1: solve QP subproblem equ (7) for primal and dual direction
        qp_solution = FunObj.qp_solver('h', qp_H_j, 'g', qp_g_j, 'a', qp_A_j, ...
            'lba', qp_lba_j, 'uba', qp_uba_j,...
            'x0', zeros(NLP.xDim + NLP.hDim + NLP.gDim, 1));
        dx = full(qp_solution.x(1 : NLP.xDim, 1));
        dlambda = full(qp_solution.x(NLP.xDim + 1 : NLP.xDim + NLP.hDim, 1)) - bar_lambda_k;
        dmu = full(qp_solution.x(NLP.xDim + NLP.hDim + 1 : end, 1)) - bar_mu_k;
        
        qp_status = FunObj.qp_solver.stats.return_status;
        disp(['inner iter j = ', num2str(j), ' --> qpoases return status: ', qp_status])
        if strcmp(qp_status, 'Successful return.')
            % step 2: check whether we can obtain a new iterate with full step along the direction
            check_lambda = (min(bar_lambda_k + dlambda) >= Option.bar_lambdaMin) && ...
                (max(bar_lambda_k + dlambda) <= Option.bar_lambdaMax) ;
            check_mu = (min(bar_mu_k + dmu) >= 0) && (max(bar_mu_k + dmu) <= Option.bar_muMax) ;
            rho_est = FunObj.rho(hat_x_j + dx, bar_lambda_k + dlambda, bar_mu_k + dmu);
            rho_est = full(rho_est);
            equ8_is_satisfied = (check_lambda) && (check_mu) && (rho_est <= r_k);
            if (strcmp(HessianType, 'LagrangianHessian')) && (equ8_is_satisfied)
                % set new outer iterate (sSQP iteration) and then break this while loop
                has_find_new_outer_iterate = true;
                hat_x_kNext_est = hat_x_j + dx;
                bar_lambda_kNext_est = bar_lambda_k + dlambda;
                bar_mu_kNext_est = bar_mu_k + dmu;
                r_kNext_est = rho_est;
                epsilon_kNext_est = epsilon_k;
                sigma_kNext_est = Option.q * r_k;
                iterateType_kNext_est = 'sSQP_iterate';
                break
            end
            % step 3: check whether we can obtain a primal direction which is descent for augmented Lagrangian
            AugLx_eta = FunObj.AugLx_eta(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k, dx);
            AugLx_eta = full(AugLx_eta);
            norm_dx = norm(dx);
            equ11_is_satisfied = (AugLx_eta <= (-Option.gamma*norm_dx^2));
            if equ11_is_satisfied
                % break this while loop and perform step 5
                break
            end
            % step 4: modify Hessian
            if omega >= Option.omegaMax
                % break this while loop and terminate the overall rountie
                modifiedHessianFlag = false;
                break
            end
            modify_Hessian = [omega*eye(NLP.xDim),    zeros(NLP.xDim, NLP.hDim + NLP.gDim);...
                             zeros(NLP.hDim + NLP.gDim, NLP.xDim + NLP.hDim + NLP.gDim)];
            qp_H_j = qp_H_j + modify_Hessian;
            HessianType = 'ModifiedHessian';
            omega = 10 * omega;
        else
            % step 4: modify Hessian
            if omega >= Option.omegaMax
                % break this while loop and terminate the overall rountie
                modifiedHessianFlag = false;
                break
            end
            modify_Hessian = [omega*eye(NLP.xDim), zeros(NLP.xDim, NLP.hDim + NLP.gDim);...
                             zeros(NLP.hDim + NLP.gDim, NLP.xDim + NLP.hDim + NLP.gDim)];
            qp_H_j = qp_H_j + modify_Hessian;
            HessianType = 'ModifiedHessian';
            omega = 10 * omega;
        end
    end
    
    %% Algorithm 1 (step 5 - step 7 if loop)
    % perform AugL iteration after step 3, i.e., the case that step 2 does not satisfy
    if (~has_find_new_outer_iterate) && modifiedHessianFlag
        % step 5: line search
        for i = 0 : Option.maxLineSearchIterNum
            % check failure flag
            if i == Option.maxLineSearchIterNum
                % break line search 'for loop'
                lineSearchFlag = false;
                break
            end
            % Amijo rules
            AugL_j = FunObj.AugL(hat_x_j, bar_lambda_k, bar_mu_k, sigma_k);
            AugL_j = full(AugL_j);
            hat_x_jNext = hat_x_j + (Option.tau)^i*dx;
            AugL_jNext = FunObj.AugL(hat_x_jNext, bar_lambda_k, bar_mu_k, sigma_k);
            AugL_jNext = full(AugL_jNext);
            if AugL_jNext <= AugL_j + Option.epsilon * (Option.tau)^i * AugLx_eta
                % terminate line search 'for loop'
                break
            end
        end
        
        if lineSearchFlag
            % step 6 check stationary of AugL
            AugLx_jNext = FunObj.AugLx(hat_x_jNext, bar_lambda_k, bar_mu_k, sigma_k);
            AugLx_jNext = full(AugLx_jNext);
            if norm(AugLx_jNext) <= epsilon_k
                % step 7: set new outer iterate (AugL iteration)
                has_find_new_outer_iterate = true;
                epsilon_kNext_est = Option.xita * epsilon_k;
                hat_x_kNext_est = hat_x_jNext;
                % dual variable
                h_jNext = FunObj.h(hat_x_jNext);
                h_jNext = full(h_jNext);
                bar_lambda_kNext_est = bar_lambda_k + 1/sigma_k*h_jNext;
                if (min(bar_lambda_kNext_est) >= Option.bar_lambdaMin) && (max(bar_lambda_kNext_est) <= Option.bar_lambdaMax)
                    % bar_lambda_kNext_est satisfies dual bounds
                else
                    bar_lambda_kNext_est = min([max([Option.bar_lambdaMin*ones(NLP.hDim, 1), bar_lambda_kNext_est], [], 2),...
                        Option.bar_lambdaMax *ones(NLP.hDim, 1)], [], 2);
                end
                g_jNext = FunObj.g(hat_x_jNext);
                g_jNext = full(g_jNext);
                bar_mu_kNext_est = max([zeros(NLP.gDim, 1), bar_mu_k + 1/sigma_k*g_jNext], [], 2);
                if (min(bar_mu_kNext_est) >= 0) && (max(bar_mu_kNext_est) <= Option.bar_muMax)
                    % bar_mu_kNext_est satisfies dual bounds
                else
                    bar_mu_kNext_est = min([max([zeros(NLP.gDim, 1), bar_mu_kNext_est], [], 2),...
                        Option.bar_muMax *ones(NLP.gDim, 1)], [], 2);
                end
                % natural residual tol and dual stabilization variable
                rho_kNext_est = FunObj.rho(hat_x_kNext_est, bar_lambda_kNext_est, bar_mu_kNext_est);
                rho_kNext_est = full(rho_kNext_est);
                if rho_kNext_est <= r_k
                    r_kNext_est = Option.q * r_k;
                    sigma_kNext_est = rho_kNext_est;
                else
                    r_kNext_est = r_k;
                    psi_kNext_est = FunObj.psi(hat_x_kNext_est, bar_mu_kNext_est);
                    psi_kNext_est = full(psi_kNext_est);
                    psi_k = FunObj.psi(hat_x_k, bar_mu_k);
                    psi_k = full(psi_k);
                    if psi_kNext_est <= Option.delta*psi_k
                        sigma_kNext_est = sigma_k;
                    else
                        sigma_kNext_est = Option.kappa*sigma_k;
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

