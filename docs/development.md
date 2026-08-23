# Development and validation

## Clean developer build

Use an out-of-source directory. In-source builds are rejected by the project.

Presets cover the usual configurations:

```sh
cmake --preset gcc-14
cmake --build --preset gcc-14
ctest --preset gcc-14
```

`clang-18`, `debug`, `asan` and `docs` are also available. Configuring by hand
works too:

```sh
cmake -S . -B build \
  -DINTERPOLATION_BUILD_EXAMPLES=ON \
  -DINTERPOLATION_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The top-level defaults enable examples and tests. When Interpolation is
consumed with `add_subdirectory` or `FetchContent`, developer targets are not
added. Pass `-DINTERPOLATION_WARNINGS_AS_ERRORS=ON` to reproduce the CI
warning policy.

## API documentation

```sh
cmake -S . -B build -DINTERPOLATION_BUILD_DOCS=ON
cmake --build build --target InterpolationDocs
```

Doxygen is required only when `INTERPOLATION_BUILD_DOCS` is enabled. Graphviz is
not required. The generated entry point is `build/docs/html/index.html`.
Doxygen warnings are treated as errors, and generated files stay in the build
tree.

## Testing conventions

- Add deterministic known-answer data for every new regression.
- Use a clearly nonuniform grid when interval indexing matters.
- Exercise real and complex ordinates where the public concepts permit both.
- Scale floating-point comparisons by the value magnitude and
  `1000 * std::numeric_limits<T>::epsilon()`.
- Keep independent reference calculations separate from production internals.
- Do not regenerate expected numerical output merely because behavior changed;
  explain and review every oracle change.
- Randomised checks take an explicit seed and print it, so a failure can be
  replayed. Override with `INTERPOLATION_TEST_SEED` to run them on different
  data.

## Benchmarks

The algorithmic changes made here have a harness that measures them against
the implementations they replaced, rather than resting on complexity
arguments. It covers barycentric `Lagrange`, the tridiagonal spline solve, and
the shared-factorisation path for a nodal derivative over many lines on one
grid:

```sh
cmake -S . -B build -DINTERPOLATION_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target InterpolationBenchmarks
./build/benchmarks/InterpolationBenchmarks
```

It is off by default and never run in CI, since timings on a shared runner
are noise. It reports the best of several repetitions and the largest relative
difference between the two implementations, so a speedup that changed the
answers is visible rather than hidden.

## Continuous integration

CI builds the gcc-13, gcc-14 and clang-18 matrix across Release and Debug with
warnings as errors, runs an ASan and UBSan job, checks formatting with
clang-format 18, builds the documentation, and installs the package to verify
that a separate consumer project can find it through `find_package`.

GCC 13 is in the matrix because it is the compiler on the deployment target
for the codes the consuming libraries serve. Nothing needs fixing for it; the
leg exists so that it stays that way by construction rather than by
inspection. The one decision that would break it is deducing `this`, which
GCC 13 lacks and which no header uses.

Reproduce any of these locally with the matching preset before pushing.

Finish each phase with `git diff --check` and review `git status --short` to
ensure generated files and unrelated repositories have not changed.
