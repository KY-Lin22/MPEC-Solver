function showInfo(self)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here


%% problem information
Dim = self.MPEC.Dim;

disp('*------------------------ MPEC Problem Information --------------------------------*')
% number of variable
disp(['Number of Total Variables: ........................', num2str(Dim.Y)])
disp(['Number of Primal Variable: ........................', num2str(Dim.Z)])
disp(['- x(optimal variable): ............................', num2str(Dim.x)])
disp(['- p(algebraic variable): ..........................', num2str(Dim.p)])
disp(['- w(auxiliary variable): ..........................', num2str(Dim.w)])

disp(['Number of Dual Variable: ..........................', num2str(Dim.LAMBDA)])
disp(['- sigma(inequality): ..............................', num2str(Dim.sigma)])
disp(['- eta(equality): ..................................', num2str(Dim.eta)])
disp(['- gamma(inequality for equilibrium constraint): ...', num2str(Dim.gamma)])
disp('')

%% solver option information
Option = self.Option;

disp('*------------------------ NIPMPEC Solver Option Information -----------------------*')
disp('1. Basic Options')
disp(['- maxIterNum: .............................', num2str(Option.maxIterNum)])
disp(['- KKT Error Tolerance(Total): .............', num2str(Option.Tolerance.KKT_Error_Total)])
disp(['                     (Feasibility): .......', num2str(Option.Tolerance.KKT_Error_Feasibility)])
disp(['                     (Stationarity): ......', num2str(Option.Tolerance.KKT_Error_Stationarity)])
disp('2. Options for Function and Jacobian Evaluation')
disp(['- Employ Sparse Pattern:...................', mat2str(Option.employSparsePattern)])
disp(['- Hessian Approximation Method: ...........', Option.HessianApproximation])
disp(['- Singularity Regular Parameter(J): .......', num2str(Option.RegularParam.nu_J)])
disp(['                               (G): .......', num2str(Option.RegularParam.nu_G)])
disp(['                               (H): .......', num2str(Option.RegularParam.nu_H)])
disp('3. Options for Search Direction Evaluation')
if Option.employSparsePattern
    disp(['- Solve Linear System Method: .............', 'mldivide_sparse'])
else
    disp(['- Solve Linear System Method: .............', Option.linearSystemSolver])
end
disp('4. Options for Line Search')
disp(['- Employ Second Order Correction: .........', mat2str(Option.employSecondOrderCorrection)]);
disp(['- stepSize_Min: ...........................', num2str(Option.LineSearch.stepSize_Min)])
disp('5. Options for Feasibility Restoration Phase')
disp(['- Employ Feasibility Restoration Phase: ...', mat2str(Option.employFeasibilityRestorationPhase)])
disp(['- maxIterNum: .............................', num2str(Option.FRP.maxIterNum)]);
disp(['- stepSize_Min: ...........................', num2str(Option.FRP.stepSize_Min)]);
disp('6. Options for Perturbed Parameter')
disp(['- sInit: ..................................', num2str(Option.sInit)])
disp(['- sEnd: ...................................', num2str(Option.sEnd)])
disp(['- zInit: ..................................', num2str(Option.zInit)])
disp(['- zEnd: ...................................', num2str(Option.zEnd)])
disp('*----------------------------------------------------------------------------------*')
end

