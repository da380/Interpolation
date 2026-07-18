# Interp

C++ header-only template library for interpolation. Linear algebra is provided
by Eigen3.

## Cubic spline boundary conditions

`CubicSplineBC::Free` is the natural boundary condition: the endpoint second
derivative is zero. `CubicSplineBC::Clamped` instead specifies the endpoint
first derivative.

The spline coefficient solve uses `Eigen::SimplicialLDLT` on a symmetric
positive-definite system stored through its lower triangle. Couplings adjacent
to Free endpoints are eliminated because those endpoint second derivatives are
known to be zero.
