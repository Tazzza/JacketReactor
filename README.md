# JacketReactor

A solution to the tubular reactor model as described in [Keßler et al](https://www.sciencedirect.com/science/article/pii/S0098135416303738). This is an example of Multi-Objective Optimal Control (MOOC) problem, where the jacket temperature is the decision variable and is encoded with different solution representations. 

This code uses the Plant Propagation Algorithm [Fresa](https://www.ucl.ac.uk/~ucecesf/fresa.html) to optimise the problem. Additionally, a multiple simultaneous solution representation is explored as explored by [Eric S. Fraga](https://arxiv.org/pdf/2106.05096).

## Solution Representations

This code currently uses two solution representations, a piecewise linear approach and a piecewise polynomial approach.

### Piecewise Linear

The piecewise linear framework is simply a uniform discretisation of the length axis, controlled by a `global const N`.

The level of discretisation can be adjusted in the `JacketReactor.jl` file.

### Piecewise Polynomial

The piecewise polynomial framework splits the length axis into sub-domains, filling each with a cubic polynomial. These temperature values and gradients are matched between sub-domains, resulting in a smooth, continuous function. 
