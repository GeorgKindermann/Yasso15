# Yasso15
Yasso15 soil carbon model
Includes also the Yasso07 Fortran code and coefficients.

The core calculaton, located in y15_subroutine.f90, is rewirtten in C++.

The linear solver, based originaly on gaussian elimination, is repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms is selectable.

To the original method, calculating matrix exponential with Taylor series, a method directly computing exp(A * t) * v using Chebyshev from Expokit was added.

A simple example shows the usage inside R.

