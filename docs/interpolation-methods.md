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
- `BoundaryCondition::NotAKnot` makes the third derivative continuous across
  the first and last interior knots, so the first two pieces are one cubic
  and the last two another. It imposes nothing false at the boundary, needs
  no extra information, and reproduces a cubic exactly. It constrains the
  whole system rather than one endpoint, so it applies at both ends and
  needs at least four nodes. It is solved in the nodal-slope formulation
  and converted back to curvatures: written directly in curvatures its
  boundary row leaves the tridiagonal band, and eliminating that entry
  produces a leading coefficient that vanishes on a uniform grid.

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
| `Quotient` | Reciprocal recurrence from the Leibniz rule | all |
| `Composition` | Faa di Bruno, via partial Bell polynomials | all |
| `DerivativeNode<K>` | Index shift onto the operand | all the operand allows |
| `PiecewisePrimitive` | Operand shifted down one order | all the operand allows |

Every node differentiates to any order. The quotient uses
\f$q^{(n)} = \left(f^{(n)} - \sum_{k<n}\binom{n}{k}q^{(k)}g^{(n-k)}\right)/g\f$,
which follows from \f$f = qg\f$; the composition uses Faa di Bruno in the
partial Bell polynomial form, built by the recurrence
\f$B_{n,k} = \sum_i \binom{n-1}{i-1} g^{(i)} B_{n-i,k-1}\f$ rather than a sum
over set partitions. Both gather the operands' derivatives once and run a
small table over them, so neither allocates.

The one remaining ceiling is `Lagrange`, which evaluates only the value and
first derivative and refuses to compile above that rather than returning a
plausible wrong number.

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

The edge condition applies on all four edges and defaults to not-a-knot, which
keeps the scheme fourth order there. The natural condition is also available
and needs only three nodes per axis, but it fixes the second derivative to zero
at the edges, which the data rarely satisfies, and the resulting
\f$O(h^2)\f$ boundary error drags global convergence down to second order.
Measured on \f$\sin x \cos y\f$, worst error over the whole square:

| grid | natural | not-a-knot |
| --- | --- | --- |
| 9 x 9 | \f$7.0\times10^{-3}\f$ | \f$5.1\times10^{-4}\f$ |
| 17 x 17 | \f$1.7\times10^{-3}\f$ | \f$3.4\times10^{-5}\f$ |
| 33 x 33 | \f$4.3\times10^{-4}\f$ | \f$2.2\times10^{-6}\f$ |
| 65 x 65 | \f$4.2\times10^{-5}\f$ | \f$5.4\times10^{-8}\f$ |

Not-a-knot converges at fourth order globally; natural does not.
