# Yasso15
Yasso15 soil carbon model
Includes also the Yasso07 Fortran code and coefficients.

The core calculaton, located in y15_subroutine.f90, is rewirtten in C++ and R.

The linear solver, based originaly on gaussian elimination, is repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms is selectable.

To the original method, calculating matrix exponential with Taylor series, a method directly computing exp(A * t) * v using Chebyshev from Expokit was added.

A simple example shows the usage inside C++ and R.


# Yasso15 for R Users

Simply download [simpleExample.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso15/master/simpleExample.r) showing how to use it and if you want a local copy of the function also download [y15_subroutine.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso15/master/y15_subroutine.r).