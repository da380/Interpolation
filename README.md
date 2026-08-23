# Interpolation

Interpolation is a C++23 header-only template library for one-dimensional
interpolation and polynomial operations. Abscissae are real floating-point
values; ordinates and polynomial coefficients may be real or complex where
noted below. The library has no external dependencies.

## Available routines

| Routine | Purpose | Minimum data |
| --- | --- | --- |
| `Linear` | Piecewise-linear interpolation and derivative | 2 nodes |
| `CubicSpline` | Smooth piecewise-cubic interpolation | 2 nodes |
| `AkimaSpline` | Local piecewise-cubic interpolation for real or complex ordinates | 3 nodes |
| `Lagrange` | Global polynomial interpolation of sampled values | 1 node |
| `LagrangeBasis` | Individual Lagrange cardinal basis functions | 1 node |
| `Polynomial` | Polynomial evaluation, calculus, and arithmetic | Defaults to zero; otherwise 1 coefficient |
| `CubicSplineSystem` | The spline system for a fixed grid, factorised once and reused across datasets | 2 nodes |
| `TridiagonalFactorization` | A tridiagonal matrix reduced once, ready to act on many right-hand sides | 1 row |

Two-dimensional interpolation on rectilinear grids:

| Routine | Purpose | Minimum data |
| --- | --- | --- |
| `Bilinear` | Tensor-product linear interpolation on a grid | 2 nodes per axis |
| `BicubicSpline` | Tensor-product cubic-spline interpolation on a grid | 2 nodes per axis |

Every one of the one-dimensional routines models the `Function1D` concept, which is what lets them be
combined arithmetically, differentiated, composed and integrated as functions
rather than only evaluated at a point.

All interpolators take ranges. An lvalue container is borrowed, so it must
outlive the interpolator and must not be reallocated; an rvalue is moved in
and owned, so the interpolator can outlive its source. Note that the splines
work out their coefficients once, at construction: changing a borrowed
container afterwards leaves the new samples paired with the old coefficients,
so rebuild rather than mutate. Abscissae must be
strictly increasing, and this is enforced at construction.

## Requirements

- A C++23 compiler. GCC 13, GCC 14 and Clang 18 are the versions CI covers
- CMake 3.24 or newer
- No external dependencies. Git access is needed only when CMake fetches
  GoogleTest to build the tests
- Doxygen 1.9 or newer only when building the API documentation

The standard is carried on the exported target, so consumers do not need to
set `CMAKE_CXX_STANDARD` themselves.

## Building

```sh
cmake --preset gcc-14
cmake --build --preset gcc-14
ctest --preset gcc-14
```

`gcc-13`, `clang-18`, `debug`, `asan` and `docs` presets are also
available. The `gcc-13` leg exists because that is the compiler on the
deployment target for the codes the consuming libraries serve.

## Use from CMake

Interpolation exports the interface target `Interpolation::Interpolation`.
It can be brought into a project with `FetchContent`:

```cmake
include(FetchContent)
FetchContent_Declare(
  Interpolation
  GIT_REPOSITORY https://github.com/da380/Interpolation.git
  GIT_TAG <reviewed-commit-or-tag>
)
FetchContent_MakeAvailable(Interpolation)

target_link_libraries(my_target PRIVATE Interpolation::Interpolation)
```

Pin a reviewed commit or release tag rather than a moving branch in consuming
projects.

## Minimal example

```cpp
#include <Interpolation/CubicSpline.hpp>
#include <vector>

std::vector<double> x{0.0, 1.0, 3.0, 4.0};
std::vector<double> y{0.0, 1.0, 0.0, 2.0};

// Two ranges select a natural spline: S'' is zero at both endpoints.
// Passing lvalues borrows them; the containers must outlive the spline.
Interpolation::CubicSpline natural{x, y};
double value = natural(2.0);            // the same as natural.Evaluate<0>(2.0)
double derivative = natural.Evaluate<1>(2.0);

// Clamped conditions specify the endpoint first derivatives.
Interpolation::CubicSpline clamped{
    x, y, Interpolation::BoundaryCondition::Clamped, -0.5, 1.25};

// Passing rvalues moves the data in, so the spline owns it and can be
// returned from a function or outlive its source. Do this last: it leaves
// x and y moved-from.
Interpolation::CubicSpline owning{std::move(x), std::move(y)};
```

Construction throws `std::invalid_argument` if the two ranges differ in
length, are too short for the method, or the abscissae are not strictly
increasing.

## The function algebra

Anything modelling `Function1D` exposes `Evaluate<N>(x)`, the `N`th derivative
at `x`, with `operator()` as `Evaluate<0>`. That one change turns
differentiation from a member function that returns a number into a node that
is itself a function, so it composes:

```cpp
#include <Interpolation/Function.hpp>

Interpolation::CubicSpline s{x, y};

auto g  = Derivative(s) * s + 2.0;   // a function, not a number
auto dg = Derivative(g);             // differentiate the product
auto S  = Primitive(s);              // antiderivative, vanishing at x.front()

double value = g(1.5);
double area  = S.Integral(0.0, 3.0);
```

Every node differentiates to any order: products through the general Leibniz
rule expanded at compile time, quotients through the recurrence that follows
from `f = qg`, and compositions through Faa di Bruno in partial Bell
polynomial form. None of them allocates. The exception is `Lagrange`, which
evaluates only the value and first derivative and refuses to compile above
that rather than returning a wrong number.

Because the operators take operands by value, something expensive or
impossible to copy — a `Piecewise` whose pieces own their samples, for
instance — enters the algebra through `Ref`, which borrows explicitly:

```cpp
auto slope = Derivative(Ref(model));   // model must outlive it
```

Nodes store their operands **by value**. That is what makes
`CubicSpline{x, y} * CubicSpline{x, y}` valid: a node holding references would
dangle as soon as it was built from a temporary, which is what every
subexpression is. Interpolators are cheap to copy, since they hold views
unless they were given rvalues. Evaluation allocates nothing, and there is a
test that counts allocations to keep it that way.

Use `<Interpolation/Linear.hpp>`, `<Interpolation/AkimaSpline.hpp>`,
`<Interpolation/Lagrange.hpp>`, or `<Interpolation/Polynomial.hpp>` to include a
single facility, and `<Interpolation/Interpolation.hpp>` for the complete
public API.

## Piecewise-continuous functions

Real data is often piecewise continuous — a layered model, say, where a
quantity jumps at an interface. The base interpolators require strictly
increasing abscissae and reject such data outright. `Piecewise` is where the
discontinuity is represented instead, as structure rather than as a
coincidence in the samples.

A doubled abscissa is one common convention for flagging a jump, so there is a
factory that reads it:

```cpp
#include <Interpolation/Piecewise.hpp>

// radius has 3.0 twice; density differs across it.
auto model = Interpolation::SplitAtRepeats(radius, density,
    [](auto r, auto d) { return Interpolation::CubicSpline{std::move(r), std::move(d)}; });

auto [below, above] = model.Limits(3.0);   // both one-sided values
double value = model(3.0);                 // right-continuous by default

auto layer = model.Piece(1);               // the outer layer, on its own
double lo = layer.Lower(), hi = layer.Upper();
double mass = Primitive(layer).Integral(lo, hi);
```

Extracting a piece is usually more useful than evaluating through the whole
object: a single layer can be handed to a quadrature or an ODE solver that
will work over exactly that interval, and it carries that interval with it.
The extracted view borrows, so it is valid while the model lives.

Three deliberate choices. The pieces **tile** their interval, since a gap would
mean the function is undefined there, which is a different thing. Continuity
at a breakpoint is **not checked** — whether the pieces agree is your business,
and enforcing it would only start an argument about tolerance. And evaluation
is **right-continuous** by default, with `Limits` giving both sides, which at a
real interface is usually what is wanted.

All pieces share a type. Mixing kinds is type erasure, a separate concern, and
`AnyFunction1D` provides it: it holds any function with matching `Real` and
`Scalar` behind one type, so `Piecewise<AnyFunction1D<double, double>>` lets
one layer be a spline and the next a polynomial, with `Piecewise` unchanged
and knowing nothing about erasure. It is the only place in the library that
erases a type: construction allocates once and evaluation is a virtual call,
so reach for it when you need mixed kinds, not by default. It shares rather
than clones on copy, which keeps it copyable even when what it holds is
move-only.

## Two dimensions

`Bilinear` and `BicubicSpline` interpolate on a rectilinear grid: two axes and
a flat, row-major value range, with element `(i, j)` at `i * size(y) + j`.

```cpp
#include <Interpolation/BicubicSpline.hpp>

std::vector<double> x{0.0, 1.0, 2.0, 3.0};
std::vector<double> y{0.0, 1.0, 2.0};
std::vector<double> v(x.size() * y.size());   // row-major

Interpolation::BicubicSpline s{x, y, v};

double value = s(1.5, 0.5);
double mixed = s.Evaluate<1, 1>(1.5, 0.5);    // d2/dx dy
```

`BicubicSpline` is a genuine tensor product, not a spline fitted line by line:
the second derivatives in each variable and the mixed fourth derivative are all
computed once at construction, so evaluation is two applications of the same
segment formula the one-dimensional spline uses, and allocates nothing.

The edge condition defaults to `BoundaryCondition::NotAKnot`, which keeps the
scheme fourth order at the boundary; it needs at least four nodes on each
axis. `BoundaryCondition::Natural` is also available and works with three,
but it forces a curvature the data usually does not have, which costs an
order near the edges and leaves the global error converging at second order.
On `sin(x)cos(y)` over a 65 by 65 grid the difference is about 770x in worst
global error.

## Examples

`examples/` holds a numbered series, meant to be read in order, each assuming
what came before. They print to stdout rather than writing data files, so
running one tells you immediately whether it did the right thing.

| | |
| --- | --- |
| `01-linear-interpolation` | Ranges in, `Evaluate<N>`, extrapolation, validation |
| `02-cubic-splines` | Smoothness, and derivatives from the same object |
| `03-boundary-conditions` | Natural, clamped and not-a-knot against exact cubic data |
| `04-akima-and-lagrange` | Overshoot near a step; the global polynomial and its basis |
| `05-polynomials` | Arithmetic, calculus, type promotion, seeded randomness |
| `06-borrowing-and-owning` | What the value category of the argument decides |
| `07-function-algebra` | `Function1D`, composing, differentiating expressions |
| `08-derivatives-and-integrals` | `Primitive`, and integrals that are exact |
| `09-piecewise-layers` | Discontinuities, extracting a layer, mixed kinds |
| `10-two-dimensions` | Grids, mixed partials, why the edge condition matters |
| `11-one-grid-many-datasets` | Factorising a grid once, and evaluating at every node |

```sh
cmake --preset gcc-14
cmake --build --preset gcc-14 --target examples
./build/gcc-14/examples/03-boundary-conditions
```

## One grid, many datasets

A cubic spline's matrix depends on the nodes alone; only the right-hand side
carries the ordinates. So a caller holding many datasets on one grid — the
components of a field, an ensemble of profiles, a line of spectral
coefficients — should factorise once and solve many times rather than
construct a spline per dataset. `CubicSplineSystem` names that, and
`CubicSpline` is built on it, so there is one assembly of the spline equations
rather than one per caller.

```cpp
#include <Interpolation/CubicSplineSystem.hpp>

const auto system = Interpolation::CubicSplineSystem{
    radius, Interpolation::BoundaryCondition::NotAKnot};   // factorised here

std::vector<std::complex<double>> curvature(radius.size());
std::vector<std::complex<double>> derivative(radius.size());

for (const auto& line : lines) {
    system.Solve(line, std::span{curvature});
    system.EvaluateAtNodes<1>(line, curvature, std::span{derivative});
}
```

`Solve` and `EvaluateAtNodes` both write into storage the caller owns and are
`const`, so the loop allocates nothing and a fixed system may be shared across
threads. The matrix is real even when the ordinates are complex, so the
ordinate type is a parameter of `Solve` rather than of the class: one system
serves both. A spline that has already been built hands its system back
through `SplineSystem()`.

The other half of the same idea is `EvaluateAtNodes`, which every
one-dimensional interpolator offers: evaluating at the nodes is the one query
that needs no search, since the segment adjoining each node is known.

```cpp
const Interpolation::CubicSpline spline{x, y};
const auto slopes = spline.NodeValues<1>();               // allocating form
spline.EvaluateAtNodes<1>(std::span{buffer});             // into your storage
```

A value at a node agrees from both sides, but a high enough derivative need
not — the third for a cubic, the first for a piecewise-linear interpolant.
Evaluation is right-continuous throughout the library, so the default reports
the limit from the segment starting at the node, and `Side::Left` asks for the
other. `Piecewise` uses the same convention across its breakpoints.

## Build, test, and document

```sh
cmake -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Examples and tests are enabled by default for a top-level build. API
documentation is opt-in:

```sh
cmake -S . -B build -DINTERPOLATION_BUILD_DOCS=ON
cmake --build build --target InterpolationDocs
```

The generated API starts at `build/docs/html/index.html` and is not tracked by
Git.

## Documentation

- [Getting started](docs/getting-started.md)
- [Interpolation methods](docs/interpolation-methods.md)
- [Development and validation](docs/development.md)
- [Rebuilding roadmap](docs/roadmap.md)

## Cubic-spline boundary conditions

`BoundaryCondition::Natural` sets the endpoint second derivative to zero.
`BoundaryCondition::Clamped` specifies the endpoint first derivative.
`BoundaryCondition::NotAKnot` makes the third derivative continuous across the
first and last interior knots, which imposes nothing false at the boundary: it
reproduces a cubic exactly and keeps the scheme fourth order right to the ends,
where Natural costs an order. It needs no extra information and at least four
nodes, and it constrains the whole system, so it must be used at both ends or
neither.

The nodal second derivatives satisfy a tridiagonal system, which is solved
directly by the Thomas algorithm. The system is strictly diagonally dominant,
so no pivoting is required. Its coefficients are real even when the ordinates
are complex, so a complex spline solves a real system against a complex
right-hand side.

Natural and Clamped are assembled in the nodal-curvature formulation and
solved for what the evaluator wants. Not-a-knot is assembled in nodal slopes
and converted afterwards: written directly in curvatures its boundary row
reaches outside the tridiagonal band, and eliminating that entry leaves a
leading coefficient of `h0^2 - h1^2`, which vanishes on a uniform grid. In
slopes the diagonal of that row is `h1 > 0` for any spacing.

Both the assembled system and the tridiagonal solver behind it are public:
`CubicSplineSystem` for the spline case, `TridiagonalFactorization` and
`SolveTridiagonal` for a caller with a banded system of their own.
