clear all
clc
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
H = 2*DM.eye(2);
A = DM.ones(1,2);
g = DM.zeros(2);
lba = 10;
qp = struct;
qp.h = H.sparsity();
qp.a = A.sparsity();
S = conic('S','qpoases',qp);
disp(S)
tic
r = S('h', H, 'g', g,...
      'a', A, 'lba', lba);
toc
x_opt = r.x;
disp(x_opt)

% methods(S)
% properties(S)
% S.print_options
% S.print_option('printLevel')
% S.stats