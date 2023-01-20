# MPEC Solver

## Introduction

This repository holds the Matlab implementation of some self-written MPEC solvers. 

**NIPMPEC_CasADi** is the implementation of our previous research to the general MPEC problem formulation. Refer to the following paper dedicated to solve the continuous time optimal control problem:  

1. Kangyu Lin and Toshiyuki Ohtsuka, "A Non-Interior-Point Continuation Method for the Optimal Control Problem with Equilibrium Constraints"，2022, (under review)  <https://arxiv.org/abs/2210.10336>, 

**Stabilized_SQP_Izmailov2015** and **Stabilized_SQP_Gill2017** are two degenerate NLP solvers which use the the augmented Lagrangian to globalize the local stabilized SQP method. These two solvers are the literature reproduction of the following two paper:

1. A.F.Izmailov et al. "Combining stabilized SQP with the augmented Lagrangian algorithm", *Computational Optimization and Applications*, 2015， <https://link.springer.com/article/10.1007/s10589-015-9744-6>
2. P.E. Gill et al. "A stabilized SQP method: superlinear convergence, ", *Mathematical Programming*, 2017, <https://link.springer.com/article/10.1007/s10107-016-1066-7>

## Requirements

- MATLAB 2019b or later
- CasADi
