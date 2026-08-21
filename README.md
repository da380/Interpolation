# Interpolation

Interpolation is a C++20 header-only template library for one-dimensional
interpolation and polynomial operations. Abscissae are real floating-point
values; ordinates and polynomial coefficients may be real or complex where
noted below. Eigen provides the linear algebra used by cubic splines.

## Available routines

| Routine | Purpose | Minimum data |
| --- | --- | --- |
| `Linear` | Piecewise-linear interpolation and derivative | 2 nodes |
| `CubicSpline` | Smooth piecewise-cubic interpolation | 2 nodes |
| `Akima` | Local piecewise-cubic interpolation for real or complex ordinates | 3 nodes |
| `Lagrange` | Global polynomial interpolation of sampled values | 1 node |
| `LagrangePolynomial` | Individual Lagrange cardinal basis functions | 1 node |
| `Polynomial1D` | Polynomial evaluation, calculus, and arithmetic | Defaults to zero; otherwise 1 coefficient |

All interpolators store iterators rather than copying their input. Keep the
sample containers alive and do not reallocate them while an interpolator is in
use. Abscissae must be strictly increasing.

## Requirements

- A C++20 compiler
- CMake 3.5 or newer
- Git access when CMake fetches Eigen and, for tests, GoogleTest
- Doxygen 1.9 or newer only when building the API documentation

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

// Three arguments select a natural spline: S'' is zero at both endpoints.
Interpolation::CubicSpline natural{x.begin(), x.end(), y.begin()};
double value = natural(2.0);
double derivative = natural.Derivative(2.0);

// Clamped conditions specify the endpoint first derivatives.
Interpolation::CubicSpline clamped{
    x.begin(), x.end(), y.begin(), Interpolation::CubicSplineBC::Clamped,
    -0.5, 1.25};
```

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

`CubicSplineBC::Free` is the natural boundary condition: the endpoint second
derivative is zero. `CubicSplineBC::Clamped` specifies the endpoint first
derivative.

The coefficient solve uses `Eigen::SimplicialLDLT` on a symmetric
positive-definite system stored through its lower triangle. Couplings adjacent
to Free endpoints are eliminated because those endpoint second derivatives are
known to be zero.
