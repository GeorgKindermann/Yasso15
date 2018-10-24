# Yasso15
Yasso15 soil carbon model

The core calculaton, located in y15_subroutine.f90, will be rewirtten in C++.

The linear solver, based originaly on gaussian elimination, will be repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms will be selectable.

Added an example showing the usage inside R.
