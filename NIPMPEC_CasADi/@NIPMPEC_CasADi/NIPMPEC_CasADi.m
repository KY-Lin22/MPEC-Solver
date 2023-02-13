classdef NIPMPEC_CasADi < handle
    %implementation of our previous research to the general MPEC problem
    %formulation:
    % min L(x,p)
    %s.t. G(x,p) >= 0
    %     C(x,p) = 0
    %     p \in SOL([l, u], K(x,p))
    %
    %Refer to the following paper dedicated to solve the continuous time optimal 
    %control problem:Kangyu Lin and Toshiyuki Ohtsuka, "A Non-Interior-Point 
    %Continuation Method for the Optimal Control Problem with Equilibrium Constraints"
    %   Detailed explanation goes here
    
    properties
        MPEC % struct, MPEC problem dimmension, symbolic representation, and data record, 
             % with field 'Dim' (variable dimension),
             %           'x'(optimal variable), 
             %           'p'(algebraic variable),
             %           'w'(auxilary variable for function K)
             %           'sigma'(dual variable for G), 
             %           'eta'(dual variable for C), 
             %           'gamma'(dual variable for PHI),
             %           'l', 'u'(lower and upper bounds for p), 
             %           'K'(function used to define BVI for p), 
             %           's'(perturbed parameter for Scholtes reformulation),
             %           'z'(perturbed parameter for FB function)
             %           'L'(cost function),
             %           'G'(inequality constraint), 
             %           'C'(equality constraint), 
             %           'PHI'(Scholtes reformulation inequalities PHI >= 0)
             %           'ZRef'(reference primal variable in FRP cost function), 
             %           'ZWeight'(weight matrix in FRP cost function),
             %           'FRP_L'(FRP cost function),
             %           
        Option % struct, solver option
        FunObj  % struct, CasADi function object 
                % with field 'K',
                %            'L', 'G', 'C', 'PHI',
                %            'L_grad', 'G_grad', 'C_grad', 'PHI_grad',
                %            'LAG_hessian', 'L_hessian',
                %            'FB_G', 'FB_G_grad', 'FB_PHI','FB_PHI_grad',
                %            'FRP_L', 'FRP_L_grad', 'FRP_L_hessian'
    end
    %% Constructor Method for NIPMPEC_CasADi  
    methods
        function self = NIPMPEC_CasADi(MPEC)
            %NIPMPEC_CasADi: Construct an instance of this class
            % Syntax:
            %          self = NIPMPEC_CasADi(MPEC)
            % Argument:
            %          MPEC: struct, containing MPEC problem formulation, 
            %                with field: 
            %                'x' -- SX symbolic variable, optimal variable
            %                'p' -- SX symbolic variable, algebraic variable
            %                'L' -- SX symbolic variable, cost function L(x,p)
            %                'G' -- SX symbolic variable, inequality constraint G(x,p) >= 0
            %                'C' -- SX symbolic variable, equality constraint C(x,p) = 0
            %                'K' -- SX symbolic variable, function used to define Box constraint Variational Inequality(BVI) for p
            %                'l' -- double, lower bounds for p
            %                'u' -- double, upper bounds for p
            %
            % Output:
            %          self: instance of this class
            
            % import CasADi to workspace
            addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
            import casadi.* 
            %% check and poblish input
            MPEC = self.checkInput(MPEC);
            
            %% initialize properties: MPEC       
            % symbolic representation of problem variable and function
            Dim = struct('x', size(MPEC.x, 1), 'p', size(MPEC.p, 1), 'w', size(MPEC.p, 1));         
            self.MPEC = struct('Dim', Dim,...
                'x', MPEC.x , 'p', MPEC.p, 'w', SX.sym('w', Dim.w, 1),...
                'l', MPEC.l, 'u', MPEC.u, 'K', MPEC.K,...
                's', SX.sym('s', 1, 1), 'z', SX.sym('z', 1, 1));                          
            self.MPEC.Dim.Z = self.MPEC.Dim.x + self.MPEC.Dim.p + self.MPEC.Dim.w; 
            
            % problem reformulation: L G C
            pWeight = 0.001 * eye(self.MPEC.Dim.p);
            wWeight = 0.001 * eye(self.MPEC.Dim.w);           
            self.MPEC.L = MPEC.L ...
                + 0.5 * self.MPEC.p' * pWeight * self.MPEC.p ...
                + 0.5 * self.MPEC.w' * wWeight * self.MPEC.w;           
            self.MPEC.G = MPEC.G;
            self.MPEC.C = [MPEC.C;...
                self.MPEC.w - self.MPEC.K];     
            
            % problem reformulation: PHI
            NCP_num = 0;
            for i = 1 : self.MPEC.Dim.p
                if (self.MPEC.l(i) == 0) && (self.MPEC.u(i) == Inf)
                    NCP_num = NCP_num + 1;
                end
            end
            BVI_num = self.MPEC.Dim.p - NCP_num;            
            PHI = SX.sym('PHI', 3 * NCP_num + 4 * BVI_num, 1);
            PHI_Counter = 0;
            for i = 1 : self.MPEC.Dim.p
                if (self.MPEC.l(i) == 0) && (self.MPEC.u(i) == Inf)
                    % nonlinear complementary problem
                    PHI(PHI_Counter + 1 : PHI_Counter + 3, 1) = ...
                        [self.MPEC.p(i);...
                        self.MPEC.w(i);...
                        self.MPEC.s - self.MPEC.p(i) * self.MPEC.w(i)];                    
                    PHI_Counter = PHI_Counter + 3;
                else
                    % box constraint variation inequality
                    PHI(PHI_Counter + 1 : PHI_Counter + 4, 1) = ...
                        [self.MPEC.p(i) - self.MPEC.l(i);...
                        self.MPEC.u(i) - self.MPEC.p(i);...
                        self.MPEC.s - (self.MPEC.p(i) - self.MPEC.l(i)) * self.MPEC.w(i);...
                        self.MPEC.s + (self.MPEC.u(i) - self.MPEC.p(i)) * self.MPEC.w(i)];              
                    PHI_Counter = PHI_Counter + 4;
                end
            end     
            self.MPEC.PHI = PHI;
            
            % dual variable
            self.MPEC.Dim.sigma = size(self.MPEC.G, 1);
            self.MPEC.Dim.eta = size(self.MPEC.C, 1);
            self.MPEC.Dim.gamma = size(self.MPEC.PHI, 1);           
            
            self.MPEC.sigma = SX.sym('sigma', self.MPEC.Dim.sigma, 1);
            self.MPEC.eta = SX.sym('eta', self.MPEC.Dim.eta, 1);
            self.MPEC.gamma = SX.sym('gamma', self.MPEC.Dim.gamma, 1);
            
            self.MPEC.Dim.LAMBDA = self.MPEC.Dim.sigma + self.MPEC.Dim.eta + self.MPEC.Dim.gamma;
            self.MPEC.Dim.Node = cumsum([self.MPEC.Dim.sigma, self.MPEC.Dim.eta, self.MPEC.Dim.gamma,...
                self.MPEC.Dim.x, self.MPEC.Dim.p, self.MPEC.Dim.w]);
            self.MPEC.Dim.Y = self.MPEC.Dim.Z + self.MPEC.Dim.LAMBDA;
            
            % feasibility restoration phase 
            self.MPEC.ZRef = SX.sym('ZRef', self.MPEC.Dim.Z, 1);
            self.MPEC.ZWeight = SX.sym('ZWeight', self.MPEC.Dim.Z, 1);
            Z = [self.MPEC.x; self.MPEC.p; self.MPEC.w];
            self.MPEC.FRP_L = 0.5 * (Z - self.MPEC.ZRef)' * diag(self.MPEC.ZWeight) * (Z - self.MPEC.ZRef);
            
            %% initialize properties: Option
            self.Option = self.createOption();                        
            
            %% initialize properties: FunObj
            self.FunObj = self.createFunObj();
            
        end
        
    end 
    
    %% Other Methods for NIPMPEC_CasADi   
    methods
        %%
        MPEC = checkInput(self, MPEC)       
        
        Option = createOption(self)
        
        FunObj = createFunObj(self)
        
        showInfo(self) 
        
        generateInitialGuess(self)
        
        [solution, Info] = solveMPEC(self, VarInit)
        
        showResult(self, Info) 
        
        %% Methods in solveMPEC method
       % function, jacobian and Hessian evaluation
        Fun = FunctionEvaluation(self, Var, s, z, mode, FRP)          
        
        Jac = JacobianEvaluation(self, Var, s, mode, FRP)
        
        Hessian = HessianEvaluation(self, Var, Jac, s, mode, FRP)
        
        % KKT evaluation
        [KKT_Residual, KKT_Error] = computeKKT_Residual_Error(self, Var, Fun, Jac)
        
        KKT_Matrix = computeKKT_Matrix(self, Fun, Jac, Hessian)
        
        % evaluate search direction
        [dY, Info] = SearchDirection(self, KKT_Residual, KKT_Matrix) 
        
        % Line Search  
        [Var_LS, Info] = LineSearch_Merit(self, Var, Fun, Jac, beta, s, z,...
            KKT_Residual, KKT_Matrix, dY_k, mode, FRP)
        
        % Second Order Correction
        Var_SOC = SecondOrderCorrection(self, Var, Fun, KKT_Residual, KKT_Matrix, Fun_full)
        
        % Feasibility Restoration Phase
        [Var_FRP, Info] = FeasibilityRestorationPhase_MinDeviation(self, Var_Ref, Fun_Ref, Jac_Ref, s, z)
        
         % compute perturbed parameter s and z
        [s_k, z_k, Fun_k] = computePerturedParam(self, Var_k, Fun_k, s, z);
        
        % examine solution
        Info = solutionExaminer(self, solution, Record)
    end
       
end

