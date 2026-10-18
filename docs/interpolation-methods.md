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
- `BoundaryCondition::Periodic` joins the final node to the first, with the
  value and first two derivatives continuous across the seam. The final node
  closes the period, and its ordinate must equal the first exactly. Like
  not-a-knot it applies at both ends and needs at least four nodes, counting
  the closing one; it has its own section below.

The three-argument constructor is Free/Free, commonly called a natural spline.
Mixed Free/Clamped combinations are supported. The implementation assembles the
tridiagonal system for the nodal second derivatives and solves it with the
Thomas algorithm, without pivoting, which is safe because the system is
strictly diagonally dominant.

### The system is separable from the data

That matrix depends on the nodes alone; only the right-hand side carries the
ordinates. `CubicSplineSystem` is that matrix, assembled and eliminated once,
with `Solve` applying it to any number of ordinate sets. `CubicSpline` holds
one and is solved by it, so the spline equations are assembled in exactly one
place, and `BicubicSpline` builds one system per axis rather than one per line.

The separation is worth naming because "one grid, many datasets" is a common
shape rather than a special case: the components of a vector field, an
ensemble of profiles, a line of spectral coefficients at every point of a
sphere. A caller in that position pays for the factorisation once. `Solve`
allocates nothing and is `const`, so a fixed system can be shared across
threads, and because the matrix is real even for complex ordinates the
ordinate type is a parameter of `Solve` rather than of the class.

`TridiagonalFactorization` sits underneath and is public in its own right, for
a caller who assembles a different banded system; `SolveTridiagonal` remains
for the genuinely one-shot case, where keeping the pivots would be waste.

### Evaluating at the nodes

Every one-dimensional interpolator offers `EvaluateAtNodes<N>` and
`NodeValues<N>`. Evaluating at a node is the one query that needs no search,
since the segment adjoining each node is known, so the nodal sweep costs a
single pass where the general path would cost `n log n` lookups. That is the
other half of what a differentiation operator on a fixed grid wants, and
`CubicSplineSystem` offers it too, taking the ordinates and curvatures rather
than holding them.

A value at a node agrees from both sides, as do the first and second
derivatives of a cubic spline; the third jumps, as does the first derivative of
a piecewise-linear interpolant. `Side` chooses which limit is reported, and
defaults to the right-hand one, so the nodal sweep agrees with an ordinary
query at the same abscissa. This is the same convention `Piecewise` uses across
its breakpoints, and it is deliberately the only one in the library.

### The periodic condition

Data on a circle has no ends, so any end condition invents something. The
periodic spline instead treats every node as interior, with indices taken
modulo the \f$m\f$ distinct nodes, \f$h_{-1} = h_{m-1}\f$ and
\f$\delta_{-1} = \delta_{m-1}\f$:

\f[
\frac{h_{i-1}}{6} M_{i-1} + \frac{h_{i-1} + h_i}{3} M_i + \frac{h_i}{6} M_{i+1}
  = \delta_i - \delta_{i-1}, \qquad i = 0, \ldots, m - 1.
\f]

That system is cyclic: tridiagonal but for one entry in each far corner. It is
solved by Sherman-Morrison. Folding the corners into the first and last
diagonal entries leaves a tridiagonal matrix that is still strictly diagonally
dominant; it is factorised once, the one auxiliary solution the rank-one
correction needs is kept with it, and each solve is then one sweep and one
correction, `const` and free of allocation as before. The correction divides
by the ratio of the two matrices' determinants, which is positive since both
are positive definite, so it never fails.

The caller passes the closing node, one period after the first, with a copy
of the first ordinate. That is how a periodic grid is usually held already,
and it is how scipy reads periodic data too. A closing ordinate that differs
by rounding is refused rather than quietly overwritten, since it is a sign the
value was computed separately.

Queries are reduced into \f$[x_0, x_n)\f$, so the spline repeats. Everything
else that answers a query follows: `Segment` reports the segment of the
reduced query; the nodal sweep wraps at the seam, so that the right-hand limit
at the closing node and the left-hand limit at the first are taken from the
pieces on the far side; and `Primitive` adds the integral over one period for
each period crossed, since the antiderivative of a function with nonzero mean
grows rather than repeating.

On equal and on stretched grids, through samples of \f$\sin\f$ and
\f$\cos\f$, the periodic spline converges at fourth order: measured rates
4.00 to 4.10.

### Exterior evaluation

Use a cubic spline for smooth material profiles or other data where continuous
first and second derivatives matter. Exterior evaluation continues the first
or final cubic piece, except under the periodic condition, and should not be
interpreted as a physically constrained extrapolation model.

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

The short form takes one condition for all four edges and defaults to
not-a-knot, which keeps the scheme fourth order there. The natural condition is
also available and needs only two nodes per axis, but it fixes the second derivative to zero
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

### Conditions per axis

Each axis can instead take its own pair of conditions, through
`AxisConditions`, with the slopes for any clamped edge in `EdgeSlopes`: for a
clamped edge of the first axis, \f$\partial_x S\f$ there at every node of the
second, and the reverse. The assembly keeps its shape, with the end data each
solve now needs:

- \f$M^x\f$, along the first axis, each line with its own two end slopes;
- \f$M^y\f$, along the second axis, under that axis's conditions;
- \f$M^{xy}\f$, along the first axis applied to \f$M^y\f$. On a clamped edge
  its end data is \f$\partial_x \partial_y^2 S\f$ there, which is the
  second-axis curvature of that edge's slopes, found by applying the
  second-axis system to the slope array.

That last step is what keeps the result the tensor product of the two
conditions. Along a clamped edge, the derivative across it is the spline of the
given slopes under the other axis's conditions, between the nodes as well as
at them, and the whole surface agrees to rounding with the spline of splines:
a one-dimensional spline along every line of the first axis, read at the
query, then one through those readings along the second.

One axis may be clamped, not both. Clamping both would need the
cross-derivative \f$\partial_x \partial_y S\f$ at the four corners as well,
which there is no natural way to supply, so the combination is refused rather
than given invented corner data.

A periodic axis follows the one-dimensional convention. Its final node closes
the period; the values along it, and the slopes on any clamped edge running
along it, must close exactly; and queries along it are reduced. The closure is
checked before any solve, so a fault is reported against `BicubicSpline`
rather than surfacing from a line solve inside the assembly; every array the
assembly derives from the data then closes by construction. One periodic
condition on all four edges gives a spline on the torus.

Clamped across the poles of \f$[0, \pi] \times [0, 2\pi)\f$ with exact slopes,
and periodic in the second variable, the scheme converges at fourth order:
measured rates 3.96 to 4.19.
