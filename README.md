# MPEC Solver

## Introduction

This repository holds the Matlab implementation of some self-written MPEC solvers. 

**NIPMPEC_CasADi** is an general-MPEC-implementation of my previous research dedicated to the continuous time optimal control problem:  

1. Kangyu Lin and Toshiyuki Ohtsuka, "A Non-Interior-Point Continuation Method for the Optimal Control Problem with Equilibrium Constraints, " <https://arxiv.org/abs/2210.10336>, 2022, (under review)

**Stabilized_SQP_Izmailov2015** and **Stabilized_SQP_Gill2017** are two degenerate NLP solvers using stabilized SQP methods with globalization scheme. These two solvers are the literature reproduction of the following two paper:

1. A.F.Izmailov et al. "Combining stabilized SQP with the augmented Lagrangian algorithm, " , <https://link.springer.com/article/10.1007/s10589-015-9744-6>, 2015
2. P.E. Gill et al. "A stabilized SQP method: superlinear convergence, ", <https://link.springer.com/article/10.1007/s10107-016-1066-7>, 2017

## Requirements

- MATLAB 2019b or later
- CasADi
