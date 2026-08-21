# Rebuilding Interpolation

Modernising this 1D header-only numerics library: zero dependencies, ranges
in, a composable function algebra, and seams cut for 2D.

Work happens on branch `refactor`. Nothing reaches `main` until the whole
sequence is proven.

A rendered copy of this plan is published at
<https://claude.ai/code/artifact/44513ac8-3970-4d9c-a6ab-c19302e2ecad>.

## Decisions

| Decision | Choice |
| --- | --- |
| Language standard | C++23, excluding `std::expected` and `std::mdspan` |
| Eigen | Removed entirely (done, phase 2) |
| Concepts | Adopt `da380/NumericConcepts`; delete the in-tree `Concepts.h` |
| Function algebra | `Function1D` concept plus free operators; no CRTP base |
| Dimensionality | 1D now, with seams cut so 2D adds API rather than breaking it |
| Namespace | Stays `Interpolation`; the README prose changes to match |
| Formatting | `.clang-format` cleaned, whole tree reformatted in one isolated commit |
| Benchmarks | Minimal in-house `std::chrono` target, opt-in, run by hand |
| Duplicate abscissae | Rejected at the base; piecewise continuity becomes its own later layer |

API breakage is acceptable: the two consuming libraries pin commit SHAs, so
they are unaffected by anything on `refactor` and can migrate deliberately.

## Where things stand

Roughly 1,400 lines of headers across `Linear`, `CubicSpline`, `Akima`,
`Lagrange` and `Polynomial1D`. The library builds clean and all 34 tests
pass, so there is a real behavioural oracle to refactor against.
Documentation and tests are in better shape than the core they describe.

Two experimental branches feed into this plan:

- `origin/developDA` holds the expression-template work: a CRTP
  `Function<Derived>` base with `Evaluate<N>(x)`, plus `Sum`, `Product` and
  `Derivative` nodes over `std::ranges` views.
- `origin/newBuildInterface` holds the `include/` layout and the move to
  NumericConcepts.

Both contain two defects to fix rather than copy. `Traits::Scalar` is
`std::invoke_result<...>` with the `_t` missing, so `Scalar` names the trait
struct rather than the type, meaning `Product`'s `Scalar{2}` cannot ever have
compiled. And `Sum` and `Product` store operands as `Derived0&`, so `f + g`
on temporaries dangles. Value-semantic nodes fix the second by construction.

## What the toolchain actually allows

Measured on the development machine rather than assumed, because it changes
the design. Both `g++-14` and `clang++-18` are installed; the default `g++`
is 13.3 and lacks deducing `this`.

| Feature | GCC 14.2 | Clang 18.1 | Verdict |
| --- | --- | --- | --- |
| Deducing `this` | yes | yes | Use it |
| `ranges::fold_left` | yes | yes | Use it |
| Ranges and views | yes | yes | Use it |
| `std::expected` | yes | no | Avoid |
| `std::mdspan` | no | no | Avoid |

Clang 18 resolves to libstdc++ 14 headers yet still cannot see
`std::expected`. The `<mdspan>` header needs GCC 15 or libc++ 18 or newer,
which is on neither toolchain.

Dropping `std::expected` costs the design nothing and keeps Clang working.
Dropping `mdspan` means the N-D seam is a small in-house extents type,
swappable for `mdspan` later without touching the public API.

## The shape of the new API

Every interpolator and expression node models a `Function1D` concept exposing
`Evaluate<N>(x)`; `operator()` is simply `Evaluate<0>`. That one change turns
differentiation from a member function into a composable node.

The highest-leverage detail is construction. Taking `viewable_range`
arguments with `views::all_t` deduction guides means an lvalue container
yields a borrowing `ref_view` and an rvalue yields an `owning_view` —
borrow-or-own from a single constructor, with no policy tag and no second
API. It also retires the dangling-iterator hazard the README currently
documents as a caveat.

```cpp
// before - three iterators, borrowed, lifetime is the caller's problem
CubicSpline s{x.begin(), x.end(), y.begin()};
double v = s(2.0);
double d = s.Derivative(2.0);

// after - ranges in; borrows lvalues, owns rvalues
CubicSpline s{x, y};                        // ref_view: borrows
CubicSpline t{std::move(x), std::move(y)};  // owning_view: owns

auto v = s(2.0);                            // same as s.Evaluate<0>(2.0)
auto d = s.Evaluate<1>(2.0);

// and the algebra falls out of the concept
auto g = Derivative(s) * s + 2.0;
auto w = g(2.0);
```

## Naming

The `value_t` / `value_type` / `x_value_type` mix goes. Abscissa precision
becomes `Real` and the ordinate type becomes `Scalar`, matching both the
`developDA` vocabulary and NumericConcepts, where `Real` and `RemoveComplex`
already mean exactly this.

This table is the full migration surface for downstream code.

| Was | Becomes | Why |
| --- | --- | --- |
| `Polynomial1D` | `Polynomial` | Dimension lives in the namespace, not the name |
| `Akima` | `AkimaSpline` | Class name matched neither its header nor its siblings |
| `LagrangePolynomial` | `LagrangeBasis` | It is the cardinal basis, not a polynomial |
| `CubicSplineBC::Free` | `BoundaryCondition::Natural` | The docs already gloss Free as "natural" |
| `x_value_type`, `value_t` | `Real` | One vocabulary across the library |
| `y_value_type` | `Scalar` | Real or complex ordinate |
| `.Derivative(x)` | `.Evaluate<1>(x)` | `Derivative(f)` is now a free function on functions |
| `<Interpolation/All>` | `<Interpolation/Interpolation.hpp>` | Extensionless headers confuse tooling |
| `.deriv()`, `.Primative()` | removed | Already deprecated forwarding aliases |

## The work, in order

The ordering is a dependency chain, not a filing scheme: each phase relies on
the one before it being green.

### Phase 0 - Safety net before anything moves

The 34 passing tests are the oracle for every later phase, so CI runs them
before a single header is touched. The students pin SHAs and are already
immune to whatever happens on `refactor`, so a tag here is courtesy rather
than insurance. Delete the committed `Testing/Temporary/` and
`output_files/*.out`, and widen `.gitignore`.

**Exit:** CI green on `refactor`; nothing pushed to `main`.

### Phase 1 - Build, packaging, CI

Raise `cmake_minimum_required` to 3.24. The current value of 3.5 sits on the
removal boundary, since CMake 4 has dropped compatibility below 3.5. Add
`VERSION`, and move to `include/Interpolation/*.hpp`.

Add `target_compile_features(INTERFACE cxx_std_23)` so consumers inherit the
standard, an `Interpolation::Interpolation` alias, and real install and
export rules so `find_package` works. NumericConcepts already has this
pattern correct and it can be lifted wholesale.

Today the include directory is the repository root, so consumers get
`tests/`, `docs/` and `examples/` on their include path. That goes. Rename
the `MY_PROJECT_BUILD_*` options, pin GoogleTest to a release tag, and add
presets. CI gets a gcc-14 and clang-18 matrix, warnings-as-errors on our own
headers, an ASan and UBSan job, a format check, and Doxygen to Pages.

**Exit:** `find_package(Interpolation)` works from an install tree; matrix
green.

### Phase 2 - Numerics core, no algebra yet

Eigen comes out. `CubicSpline` builds a sparse matrix and runs
`SimplicialLDLT` on a system that is merely tridiagonal. A Thomas solve is
about twenty lines, faster, allocation-light, and takes a real coefficient
matrix against a complex right-hand side instead of instantiating a complex
sparse matrix. `AkimaSpline.h` includes three Eigen headers and uses none of
them. That leaves the library dependency-free apart from NumericConcepts.

Delete `Concepts.h` in favour of NumericConcepts, adding only a thin local
refinement: its `RealRange` requires just `input_range`, and interpolation
needs random-access and sized ranges.

Convert all four interpolators to ranges and `Evaluate<N>`. Factor the
duplicated `upper_bound`-and-clamp into one locator. Rewrite `Lagrange` in
barycentric form; it is currently quadratic to evaluate and cubic to
differentiate, recomputing basis denominators on every call.

Fix `Polynomial`: hidden friends rather than operators in the global
namespace, `common_type_t` results instead of silently taking the left
operand's type, signed sizes, trailing-zero trimming so `Degree()` stops
lying after a subtraction, and a seedable `Random`.

**Exit:** no Eigen in the tree; ported tests green; `-Wall -Wextra` clean.

### Phase 3 - The function algebra

A `Function1D` concept plus free arithmetic operators, `Compose`,
`Derivative<K>` and `Primitive`. Nodes store operands by value, which
structurally kills the dangling-reference bug on `developDA`. Deducing
`this` supplies any shared helpers, so there is no CRTP base and no `Traits`
specialisation to write per node.

**Exit:** `Derivative(s) * s + 2.0` evaluates, allocation-free, against
analytic answers.

### Phase 4 - Second dimension

A minimal grid and extents abstraction, then tensor-product bilinear and
bicubic interpolation on rectilinear grids, reusing the 1D pieces rather
than restating them. Because the seams were cut in phases 2 and 3, this adds
API without breaking it again.

**Exit:** 2D interpolation against analytic surfaces; 1D API unchanged.

## Correctness items folded into the above

- Compiling the headers with `-Wall -Wextra` surfaces two unused variables in
  `Linear::Derivative`, dead copy-paste from `operator()`, and sign-compare
  bugs in `Polynomial1D`. The build currently sets no warning flags at all,
  so it reports none of this.
- All preconditions are `assert`, so they vanish under `NDEBUG`, including
  the check that the linear solve succeeded. New policy: throw
  `std::invalid_argument` at construction for data errors such as unsorted
  abscissae, length mismatch, or too few nodes, since construction is not the
  hot path. Evaluation stays precondition-only. No `std::expected`, so no
  compiler is excluded.
- Nothing anywhere checks that `y` is as long as `x`; only `x` has a begin
  and an end. Ranges make this checkable and cheap.
- `CubicSpline`'s default constructor leaves its iterators indeterminate, so
  evaluating such an object is undefined. It goes, or it gains a valid empty
  state.
- `Polynomial1D::Degree()` returns `size() - 1` unsigned, so an empty
  coefficient vector yields `SIZE_MAX`. This is reachable today.
- Randomised tests seed from `random_device`, so a CI failure may not
  reproduce. Seed deterministically and print the seed.

## Deferred by design

Piecewise-continuous functions are a first-class use case, not an edge case,
and they are deliberately *not* being solved by letting duplicate abscissae
through the base interpolators.

A doubled `x` is one common convention for flagging a discontinuity in
tabulated data, but it is only one of several, and honouring it implicitly
would force every evaluator to carry which-side-of-the-jump logic for the
sake of a convention the caller never stated. Instead the base requires
strictly increasing abscissae, and a later layer holds a list of base
functions together with

- a way to select the function on a particular sub-interval, and
- direct evaluation with an explicit rule for the value *at* a breakpoint,
  whether that is the left limit, the right limit, or a per-object choice.

That layer is a separate, considered piece of work after phase 4. It is
listed here so phases 2 and 3 leave room for it rather than designing it
opportunistically along the way. The relevant asymmetry to remember: a
repeated node is fatal during *setup* for barycentric `Lagrange`, which
divides by node differences, whereas `Linear` would merely see a zero-width
interval it never selects. Uniform rejection is the simpler invariant.

## Baseline

Verified 21 August 2026: configure, build and `ctest` all green, 34 of 34
tests passing, with GCC 13.3, GCC 14.2, Clang 18.1, CMake 3.28.3, and Eigen
fetched at `master`.
