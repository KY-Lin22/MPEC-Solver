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
        FunObj % struct, CasADi function object 
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
            self.FunObj = self.createFunObj();
                        
        end
        
    end
    
    
    %% Other Methods for stabilized_SQP_Izamailov2015 
    methods
        Option = createOption(self)     
        
        FunObj = createFunObj(self)
        
        [x_Opt, Info] = solveNLP(self, x_Init) 
        
        [hat_x_kNext, bar_lambda_kNext, bar_mu_kNext, r_kNext, epsilon_kNext, sigma_kNext, iterateType_kNext, innerLoopFlag] = ...
                innerLoopIteration(self, hat_x_k, bar_lambda_k, bar_mu_k, r_k, epsilon_k, sigma_k, iterateType_k)      
            
        [qp_H_j, qp_A_j, qp_g_j, qp_lba_j, qp_uba_j] = qpoases_Matrix(self, hat_x_j, bar_lambda_k, bar_mu_k, sigma_k)                
                     
    end
end

