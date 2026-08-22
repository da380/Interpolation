// 06 - Borrowing and owning
//
// Every interpolator takes ranges, and what it does with them follows from
// the value category of the argument. There is no policy flag: an lvalue is
// borrowed, an rvalue is moved in and owned.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Linear.hpp>

#include <cstdio>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

// An owning interpolator can be returned from a function. Under an
// iterator-based interface this would hand back three dangling iterators;
// here the samples come with it.
static auto
MakeOwningSpline() {
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{0.0, 1.0, 8.0, 27.0};
    return Interpolation::CubicSpline{std::move(x), std::move(y)};
}

int
main() {
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{0.0, 1.0, 8.0, 27.0};

    // Passing lvalues borrows them. The deduced view type says so.
    const Interpolation::CubicSpline borrowing{x, y};
    static_assert(
        std::is_same_v<decltype(borrowing),
                       const Interpolation::CubicSpline<
                           std::ranges::ref_view<std::vector<double>>,
                           std::ranges::ref_view<std::vector<double>>>>,
        "an lvalue container should deduce to a borrowing ref_view");

    std::printf("Borrowing: f(1.5) = %.6f\n", borrowing(1.5));
    std::printf("  The container must outlive the interpolator, and must not\n"
                "  be reallocated while it is in use.\n");

    // A borrowing interpolator reads through to the container, so a later
    // change to the samples is visible. Linear is used here on purpose: it
    // computes nothing in advance, so reading through gives the interpolant
    // of the new data.
    std::vector<double> ly{0.0, 1.0, 8.0, 27.0};
    const Interpolation::Linear line{x, ly};
    std::printf("\nBorrowing and mutating, with Linear:\n");
    std::printf("  before: f(1.5) = %.6f\n", line(1.5));
    ly[1] = 5.0;
    std::printf("  after y[1] = 5: f(1.5) = %.6f\n", line(1.5));

    // Be careful with this. An interpolator that precomputes coefficients --
    // every spline here does -- works them out once, at construction. Changing
    // the samples afterwards leaves the new values paired with the old
    // coefficients, which is neither the old interpolant nor the new one.
    // Build a fresh object instead, or hand it an rvalue so nobody else can
    // reach the data.
    std::printf("\n  Splines precompute, so mutating a borrowed container\n"
                "  after construction gives a stale, inconsistent object.\n"
                "  Rebuild instead.\n");

    // An owning interpolator has no such constraint.
    const auto owning = MakeOwningSpline();
    std::printf("\nOwning: f(1.5) = %.6f, built from vectors that no longer\n"
                "  exist, which is the case the old API could not express.\n",
                owning(1.5));

    // The same choice applies at the call site.
    const Interpolation::CubicSpline moved{std::move(x), std::move(y)};
    static_assert(
        std::is_same_v<decltype(moved),
                       const Interpolation::CubicSpline<
                           std::ranges::owning_view<std::vector<double>>,
                           std::ranges::owning_view<std::vector<double>>>>,
        "an rvalue should deduce to an owning_view");
    std::printf("\nMoved-in: f(1.5) = %.6f\n", moved(1.5));

    return 0;
}
