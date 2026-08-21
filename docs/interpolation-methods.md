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

- `BoundaryCondition::Natural` means `M=0` at that endpoint.
- `BoundaryCondition::Clamped` means `S'` at that endpoint equals the supplied
  derivative.

The three-argument constructor is Free/Free, commonly called a natural spline.
Mixed Free/Clamped combinations are supported. The implementation assembles the
tridiagonal system for the nodal second derivatives and solves it with the
Thomas algorithm, without pivoting, which is safe because the system is
strictly diagonally dominant.

Use a cubic spline for smooth material profiles or other data where continuous
first and second derivatives matter. Exterior evaluation continues the first
or final cubic piece and should not be interpreted as a physically constrained
extrapolation model.

## Akima interpolation

`AkimaSpline` creates cubic Hermite pieces whose knot derivatives are weighted from
neighboring secant slopes. It is local and can reduce oscillation around abrupt
changes compared with a global polynomial.

The implementation requires at least three nodes. At an interior knot it
selects the interval to the right; at the final knot it selects the final
interval. Exterior queries continue the first or final cubic piece. These
continuations are numerical extrapolations, not physically constrained models.
`Evaluate<N>` gives the value at `N = 0` and derivatives above that; the
pieces are cubic, so `Evaluate<3>` is piecewise constant and higher orders are
zero. Real magnitude weights allow both real and complex ordinates.

## Lagrange interpolation

For distinct nodes, `LagrangeBasis` evaluates the cardinal basis

\f[
\ell_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j}.
\f]

The basis satisfies `ell_i(x_j) = delta_ij` and its basis functions sum to one.
`Lagrange` combines these functions with sampled ordinates:

\f[
p(x)=\sum_i y_i\ell_i(x).
\f]

It is not evaluated from that product form. The implementation stores the
barycentric weights

\f[
w_j=\frac{1}{\prod_{k\ne j}(x_j-x_k)},
\f]

computed once at construction in \f$O(n^2)\f$, and evaluates the second
barycentric formula

\f[
p(x)=\left.\sum_j \frac{w_j y_j}{x-x_j}\middle/\sum_j \frac{w_j}{x-x_j}\right.
\f]

in \f$O(n)\f$ per query. Evaluating the products directly instead costs
\f$O(n^2)\f$ per value and \f$O(n^3)\f$ per derivative, recomputing the same
denominators every time.

A query that lands exactly on a node makes the quotient `0/0`, so those cases
are handled separately: the value is the sampled ordinate, and the derivative
uses the standard node formulas.

With `n` nodes, the result exactly represents polynomials of degree at most
`n-1`, apart from floating-point rounding. Evaluation is global: every sample
contributes to every query, including extrapolation.

## Polynomial operations

`Polynomial<T>` stores coefficients in ascending power order. Default
construction produces the degree-zero polynomial with coefficient `0`. It
evaluates the polynomial and its derivative with Horner-style recurrences,
evaluates the zero-constant antiderivative with `Primitive`, and computes
definite integrals with `Integrate`.

Scalar addition and subtraction modify the constant coefficient. Scalar
multiplication and division affect every coefficient. Polynomial addition and
subtraction operate coefficient-wise, and polynomial multiplication uses
coefficient convolution. The historical spelling `Primative` remains as a
deprecated forwarding alias. Cross-type construction and assignment convert
each coefficient when the source scalar type is convertible to the destination
scalar type.

## The function algebra

Every routine above models `Function1D`: it names an abscissa type `Real`, a
value type `Scalar`, and evaluates `Evaluate<N>(x)`, the `N`th derivative.

Free operators build expression nodes from those. The derivative rules are:

| Node | Derivative rule | Orders supported |
| --- | --- | --- |
| `Sum`, `Difference` | Linearity | all |
| `Product` | General Leibniz, expanded at compile time | all |
| `Quotient` | Quotient rule | value and first |
| `Composition` | Chain rule | value and first |
| `DerivativeNode<K>` | Index shift onto the operand | all the operand allows |
| `PiecewisePrimitive` | Operand shifted down one order | all the operand allows |

Where a rule is not implemented the node refuses to compile rather than
returning a plausible wrong number. Higher quotient derivatives need the
recurrence for derivatives of a reciprocal, and higher composition
derivatives need Faa di Bruno's formula.

Antidifferentiation takes two paths. A piecewise-polynomial function exposes
its nodes and a per-segment integral, and `Primitive` accumulates the
cumulative integral once when the node is built, so an interpolator carries no
antiderivative state unless one is asked for. A function with a closed-form
antiderivative, such as `Polynomial`, provides it directly. `Lagrange` supports
neither yet: it is a single global polynomial rather than a piecewise one.

## Two dimensions

`Bilinear` and `BicubicSpline` interpolate on a rectilinear grid, given two
axes and a flat, row-major value range. Both are tensor products of the
corresponding one-dimensional scheme, and both are built from the same segment
formulas the one-dimensional classes use rather than restating them.

Because a tensor product is separable, a mixed partial derivative is exactly
the two one-dimensional rules applied in turn:

\f[
\frac{\partial^{n+m} S}{\partial x^n \partial y^m}
  = B^{(m)}_y \left( B^{(n)}_x(\cdot) \right),
\f]

where \f$B_x\f$ and \f$B_y\f$ are the one-dimensional segment operators.

For `BicubicSpline` the construction solves three sets of one-dimensional
systems:

- \f$M^x\f$, the second derivative in the first variable, one solve per column;
- \f$M^y\f$, the second derivative in the second variable, one solve per row;
- \f$M^{xy}\f$, the mixed fourth derivative, obtained by applying the
  first-axis solve to \f$M^y\f$.

That last one is what makes it a true tensor product rather than two
independent one-dimensional fits. Evaluation then needs only the sixteen
corner quantities of the containing cell and allocates nothing.

The natural condition is applied on all four edges, which fixes the second
derivative to zero there. In the interior the scheme is fourth order; near the
boundary the natural condition contributes an \f$O(h^2)\f$ error, so the global
worst-case error converges at second order. Measured on
\f$\sin x \cos y\f$ over a 33 by 33 grid, the interior error is around
\f$5\times 10^{-7}\f$ against \f$1.7\times 10^{-3}\f$ for bilinear on the
same data. Not-a-knot end conditions would lift the boundary order and are the
natural next improvement.
