# Development and validation

## Clean developer build

Use an out-of-source directory. In-source builds are rejected by the project.

```sh
cmake -S . -B build \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DMY_PROJECT_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The top-level defaults enable examples and tests. When Interp is consumed with
`add_subdirectory` or `FetchContent`, developer targets are not added.

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

## Validate with Eigen 3.4.0

DSpecM1D consumes Eigen 3.4.0. Given an exact local Eigen source checkout and a
local GoogleTest source checkout, configure without editing Interp's dependency
declaration:

```sh
cmake -S . -B build-eigen-3.4 \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DMY_PROJECT_BUILD_TESTS=ON \
  -DINTERPOLATION_BUILD_DOCS=ON \
  -DBUILD_TESTING=OFF \
  -DFETCHCONTENT_SOURCE_DIR_EIGEN3=/path/to/eigen-3.4.0 \
  -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path/to/googletest
cmake --build build-eigen-3.4 --parallel
cmake --build build-eigen-3.4 --target InterpolationDocs
ctest --test-dir build-eigen-3.4 --output-on-failure
```

`BUILD_TESTING=OFF` prevents Eigen's own large test registry from being added;
Interp still enables and registers its tests through `MY_PROJECT_BUILD_TESTS`.

Finish each phase with `git diff --check` and review `git status --short` to
ensure generated files and unrelated repositories have not changed. Record the
result in `implementation_status.md`.
