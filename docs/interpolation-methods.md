# Interpolation methods

## Linear interpolation

On an interval `[x_i, x_{i+1}]`, `Linear` evaluates

\f[
L(x) = \frac{x_{i+1}-x}{h_i}y_i
     + \frac{x-x_i}{h_i}y_{i+1},
\qquad h_i=x_{i+1}-x_i.
\f]

Its derivative is the constant secant slope on the selected interval. At an
interior knot the implementation selects the interval to the right; at the
last knot it selects the final interval. Exterior queries use the first or
last segment.

Use linear interpolation when local monotonic behavior and low cost are more
important than derivative continuity.

## Cubic splines

`CubicSpline` constructs a piecewise cubic that interpolates every knot and has
continuous first and second derivatives. The solved unknowns are the knot
second derivatives `M_i = S''(x_i)`.

- `CubicSplineBC::Free` means `M=0` at that endpoint.
- `CubicSplineBC::Clamped` means `S'` at that endpoint equals the supplied
  derivative.

The three-argument constructor is Free/Free, commonly called a natural spline.
Mixed Free/Clamped combinations are supported. The implementation assembles the
lower triangle of a symmetric positive-definite system and solves it with
`Eigen::SimplicialLDLT<..., Eigen::Lower>`.

Use a cubic spline for smooth material profiles or other data where continuous
first and second derivatives matter. Exterior evaluation continues the first
or final cubic piece and should not be interpreted as a physically constrained
extrapolation model.

## Akima interpolation

`Akima` creates cubic Hermite pieces whose knot derivatives are weighted from
neighboring secant slopes. It is local and can reduce oscillation around abrupt
changes compared with a global polynomial.

The implementation requires at least three nodes. At an interior knot it
selects the interval to the right; at the final knot it selects the final
interval. Exterior queries continue the first or final cubic piece. These
continuations are numerical extrapolations, not physically constrained models.
`Derivative` is the canonical derivative method; `deriv` remains as a
deprecated compatibility alias. Real magnitude weights allow both real and
complex ordinates.

## Lagrange interpolation

For distinct nodes, `LagrangePolynomial` evaluates the cardinal basis

\f[
\ell_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j}.
\f]

The basis satisfies `ell_i(x_j) = delta_ij` and its basis functions sum to one.
`Lagrange` combines these functions with sampled ordinates:

\f[
p(x)=\sum_i y_i\ell_i(x).
\f]

With `n` nodes, the result exactly represents polynomials of degree at most
`n-1`, apart from floating-point rounding. Evaluation is global: every sample
contributes to every query, including extrapolation.

## Polynomial operations

`Polynomial1D<T>` stores coefficients in ascending power order. It evaluates
the polynomial and its derivative with Horner-style recurrences, evaluates the
zero-constant antiderivative with `Primitive`, and computes definite integrals
with `Integrate`.

Scalar addition and subtraction modify the constant coefficient. Scalar
multiplication and division affect every coefficient. Polynomial addition and
subtraction operate coefficient-wise, and polynomial multiplication uses
coefficient convolution. The historical spelling `Primative` remains as a
deprecated forwarding alias. Cross-type construction and assignment convert
each coefficient when the source scalar type is convertible to the destination
scalar type.
