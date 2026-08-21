# Getting started

## Choosing an include

Each public forwarding header exposes one facility:

| Header | Public API |
| --- | --- |
| `<Interpolation/Linear.hpp>` | `Linear` |
| `<Interpolation/CubicSpline.hpp>` | `CubicSpline`, `BoundaryCondition` |
| `<Interpolation/AkimaSpline.hpp>` | `AkimaSpline` |
| `<Interpolation/Lagrange.hpp>` | `Lagrange`, `LagrangeBasis` |
| `<Interpolation/Polynomial.hpp>` | `Polynomial` |
| `<Interpolation/Interpolation.hpp>` | All of the above and the public concepts |

Interpolation is header-only and has no external dependencies. Link the CMake
interface target `Interpolation::Interpolation` so the include directory and
the C++ standard are supplied consistently.

## Add Interp to a project

For a sibling checkout:

```cmake
add_subdirectory(path/to/Interp)
target_link_libraries(my_target PRIVATE Interpolation)
```

For a remote dependency, use the `FetchContent` example in the README and pin a
reviewed immutable commit or tag.

Interpolation installs and exports, so `find_package(Interpolation)` also
works against an install tree.

## Sample ownership and input rules

Every interpolator takes ranges, and what it does with them follows from the
value category of the argument:

- an **lvalue** container is **borrowed**. It must outlive the interpolator
  and must not be reallocated, exactly as the old iterator-based API required;
- an **rvalue** is **owned**. The data is moved in, so the interpolator can be
  returned from a function or stored beyond the lifetime of its source.

```cpp
Interpolation::CubicSpline borrowing{x, y};                       // ref_view
Interpolation::CubicSpline owning{std::move(x), std::move(y)};    // owning_view
```

Input is checked at construction, which throws `std::invalid_argument` if the
ranges differ in length, are too short for the method, or the abscissae are
not strictly increasing. Abscissae must be **strictly** increasing: a repeated
node is not an encoding for a discontinuity, and every method here divides by
a node difference somewhere. Represent a piecewise-continuous function by
composing several interpolators instead.

| Routine | Input precondition | Query behavior |
| --- | --- | --- |
| `Linear` | At least 2 nodes | Extrapolates with the first or final line segment |
| `CubicSpline` | At least 2 nodes | Extrapolates with the first or final cubic segment |
| `AkimaSpline` | At least 3 nodes | Extrapolates with the first or final cubic segment |
| `Lagrange` | At least 1 node | Evaluates the global polynomial for any real query |
| `LagrangeBasis` | At least 1 node | Evaluates the selected basis polynomial for any real query |
| `Polynomial` | Defaults to zero; explicit input needs at least 1 coefficient | Evaluates for any supported scalar argument |

## Scalar types

Abscissae must use `float`, `double`, `long double`, or a cv-qualified form of
one of those real types. Ordinates may use the same real types or
`std::complex` with a floating-point component type. The abscissa and ordinate
types must support the conversions and arithmetic required by
`InterpolationRanges`. AkimaSpline calculates its nonnegative interpolation
weights from real magnitudes, so its cubic values and derivatives support both
real and complex ordinates.

`Polynomial` supports real and complex floating-point coefficient types. Its
coefficient list is ordered from the constant term upward:

```cpp
Interpolation::Polynomial<double> p{1.0, -2.0, 3.0};
// p(x) = 1 - 2x + 3x^2
```

A default-constructed `Polynomial` is the degree-zero polynomial with the
single coefficient `0`.

Construction and assignment from a polynomial with a convertible coefficient
type convert each coefficient and replace the destination coefficient list.

## Boundary-condition constructors

The three-argument `CubicSpline` constructor selects Free conditions at both
ends. A six-argument overload applies one condition type to both endpoints,
while the full seven-argument overload selects each endpoint independently.
Derivative values passed for Free endpoints are ignored.
