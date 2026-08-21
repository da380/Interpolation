# Getting started

## Choosing an include

Each public forwarding header exposes one facility:

| Header | Public API |
| --- | --- |
| `<Interpolation/Linear.hpp>` | `Linear` |
| `<Interpolation/CubicSpline.hpp>` | `CubicSpline`, `CubicSplineBC` |
| `<Interpolation/AkimaSpline.hpp>` | `Akima` |
| `<Interpolation/Lagrange.hpp>` | `Lagrange`, `LagrangePolynomial` |
| `<Interpolation/Polynomial.hpp>` | `Polynomial1D` |
| `<Interpolation/Interpolation.hpp>` | All of the above and the public concepts |

Interp is header-only. Link the CMake interface target `Interpolation` so the
include directory, C++ standard, Eigen dependency, and transitive requirements
are supplied consistently.

## Add Interp to a project

For a sibling checkout:

```cmake
add_subdirectory(path/to/Interp)
target_link_libraries(my_target PRIVATE Interpolation)
```

For a remote dependency, use the `FetchContent` example in the README and pin a
reviewed immutable commit or tag.

Interp currently has no standalone install/export packaging. Consumers should
use `FetchContent` or `add_subdirectory`.

## Sample ownership and input rules

`Linear`, `CubicSpline`, `Akima`, `Lagrange`, and `LagrangePolynomial` retain
random-access iterators into caller-owned data. Their source containers must:

1. outlive the interpolation object;
2. remain at the same memory location; and
3. not be resized or otherwise reallocated.

The ordinate range must contain at least as many values as the abscissa range.
Abscissae are real and strictly increasing; duplicate nodes cause division by
zero in the interpolation formulas.

| Routine | Input precondition | Query behavior |
| --- | --- | --- |
| `Linear` | At least 2 nodes | Extrapolates with the first or final line segment |
| `CubicSpline` | At least 2 nodes | Extrapolates with the first or final cubic segment |
| `Akima` | At least 3 nodes | Extrapolates with the first or final cubic segment |
| `Lagrange` | At least 1 node | Evaluates the global polynomial for any real query |
| `LagrangePolynomial` | At least 1 node | Evaluates the selected basis polynomial for any real query |
| `Polynomial1D` | Defaults to zero; explicit input needs at least 1 coefficient | Evaluates for any supported scalar argument |

## Scalar types

Abscissae must use `float`, `double`, `long double`, or a cv-qualified form of
one of those real types. Ordinates may use the same real types or
`std::complex` with a floating-point component type. The abscissa and ordinate
types must support the conversions and arithmetic required by
`InterpolationIteratorPair`. Akima calculates its nonnegative interpolation
weights from real magnitudes, so its cubic values and derivatives support both
real and complex ordinates.

`Polynomial1D` supports real and complex floating-point coefficient types. Its
coefficient list is ordered from the constant term upward:

```cpp
Interpolation::Polynomial1D<double> p{1.0, -2.0, 3.0};
// p(x) = 1 - 2x + 3x^2
```

A default-constructed `Polynomial1D` is the degree-zero polynomial with the
single coefficient `0`.

Construction and assignment from a polynomial with a convertible coefficient
type convert each coefficient and replace the destination coefficient list.

## Boundary-condition constructors

The three-argument `CubicSpline` constructor selects Free conditions at both
ends. A six-argument overload applies one condition type to both endpoints,
while the full seven-argument overload selects each endpoint independently.
Derivative values passed for Free endpoints are ignored.
