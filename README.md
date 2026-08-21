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

Every one of these models the `Function1D` concept, which is what lets them be
combined arithmetically, differentiated, composed and integrated as functions
rather than only evaluated at a point.

All interpolators take ranges. An lvalue container is borrowed, so it must
outlive the interpolator and must not be reallocated; an rvalue is moved in
and owned, so the interpolator can outlive its source. Abscissae must be
strictly increasing, and this is enforced at construction.

## Requirements

- A C++23 compiler. GCC 14 and Clang 18 are the versions CI covers
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

`clang-18`, `debug`, `asan` and `docs` presets are also available.

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

Products differentiate to any order through the general Leibniz rule, expanded
at compile time. Quotients and compositions currently support the value and
first derivative and refuse to compile above that, rather than returning a
wrong number.

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

The nodal second derivatives satisfy a tridiagonal system, which is solved
directly by the Thomas algorithm. The system is strictly diagonally dominant,
so no pivoting is required. Its coefficients are real even when the ordinates
are complex, so a complex spline solves a real system against a complex
right-hand side.
