# Implementation status

## Phase 1: Doxygen API surface

Status: complete in the working tree; not committed or pushed.

Baseline: `cubic-spline-ldlt` at
`bdafbde8319f0147e5ac0fad67fd983dbec01743`.

Completed changes:

- Added the opt-in `INTERPOLATION_BUILD_DOCS` CMake option and the
  `InterpolationDocs` target.
- Added a warning-clean Doxygen 1.9.8 configuration that generates untracked
  HTML in the build directory without Graphviz.
- Documented the public concepts, interpolation classes, boundary conditions,
  iterator lifetimes, input preconditions, evaluation domains, and polynomial
  arithmetic API.
- Added `Akima::Derivative` and `Polynomial1D::Primitive`; retained `deriv` and
  `Primative` as deprecated forwarding aliases.

Validation:

- Default dependency build (Eigen 5.0.1 at
  `e0f18b799e2cc7d416fdb0aec2075dceebbc8295`): all examples compiled,
  `InterpolationDocs` completed with warnings treated as errors, and all 16
  tests passed.
- Exact Eigen 3.4.0 build (commit
  `3147391d946bb4b6c68edd901f2add6ac1f31f8c`): all examples compiled,
  `InterpolationDocs` completed with warnings treated as errors, and all 16
  tests passed.
- A disposable compile-and-run smoke check instantiated both new canonical
  template methods successfully.
- `git diff --check` passed, generated HTML remained under `/tmp`, and
  DSpecM1D remained unchanged.

The Akima limitation recorded during this phase was corrected in the combined
defect-correction sweep below.

## Phases 2 and 3: user documentation and deterministic coverage

Status: complete together in the working tree; not committed or pushed.

Completed documentation:

- Expanded the README with the feature summary, requirements, CMake
  integration, minimal cubic-spline example, build commands, and documentation
  index.
- Added getting-started, interpolation-method, and development/validation
  guides; all pages are included in the warning-strict Doxygen build.

Completed tests:

- Added deterministic nonuniform known-answer coverage for Linear, Akima's
  currently supported domain, LagrangePolynomial, Lagrange, and Polynomial1D.
- Covered values, derivatives, extrapolation where supported, polynomial
  calculus and arithmetic, real/complex data where the implementation
  supports them, public concepts, forwarding headers, and both retained
  compatibility aliases.
- Retained the existing randomized Linear and CubicSpline tests unchanged.
- Increased the registered Interp suite from 16 to 29 passing tests without
  adding fixtures or numerical-oracle files.

Combined validation:

- Default dependency build (Eigen 5.0.1 at
  `e0f18b799e2cc7d416fdb0aec2075dceebbc8295`): all examples compiled,
  `InterpolationDocs` completed with warnings treated as errors, and all 29
  tests passed.
- Exact Eigen 3.4.0 build (commit
  `3147391d946bb4b6c68edd901f2add6ac1f31f8c`): all examples compiled,
  `InterpolationDocs` completed with warnings treated as errors, and all 29
  tests passed.
- `git diff --check` passed; the original Linear/CubicSpline test sources and
  dependency declarations remained unchanged, no generated HTML entered the
  repository, and DSpecM1D remained clean.

## Combined defect-correction sweep

Status: complete in the working tree; not committed or pushed.

Completed changes:

- Clamped Akima interval selection to the first or final cubic for exterior
  queries and selected the final cubic at the final knot, removing the former
  out-of-range indexing paths.
- Defined Akima's minimum input as three nodes and corrected its three-node
  slope construction so every stored slope corresponds to one knot.
- Calculated Akima weights from real magnitudes, enabling the already-declared
  complex-ordinate API without changing the public template interface.
- Restored `Polynomial1D` copy and move assignment and made cross-type
  assignment convert coefficients into a replacement destination vector.

Added deterministic regressions for Akima's first/final knots, both exterior
regions, nonlinear three-node input, two-point rejection, and real and complex
known answers. Added compile-time and runtime checks for `Polynomial1D` copy,
move, cross-real-type, and real-to-complex assignment.

Validation:

- Fresh default dependency build (Eigen 5.0.1 at
  `e0f18b799e2cc7d416fdb0aec2075dceebbc8295`): all examples compiled, all 33
  tests passed, and warning-strict Doxygen 1.9.8 generation succeeded.
- Fresh exact Eigen 3.4.0 build (commit
  `3147391d946bb4b6c68edd901f2add6ac1f31f8c`): all examples compiled, all 33
  tests passed, and warning-strict Doxygen 1.9.8 generation succeeded.
- All 33 tests passed with libstdc++ checked iterators enabled.
- All 33 tests passed under AddressSanitizer and UndefinedBehaviorSanitizer
  without findings. Leak detection was disabled because LeakSanitizer cannot
  run under the validation environment's ptrace supervision.
- `git diff --check` passed; existing Linear/CubicSpline tests and dependency
  declarations remained unchanged, generated documentation stayed under
  `/tmp`, and DSpecM1D remained unchanged.

## Deferred issues

- Audit default-constructed `Polynomial1D` behavior and duplicate interpolation
  nodes separately; this sweep assigns neither issue new semantics.
