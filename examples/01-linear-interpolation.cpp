// 01 - Linear interpolation
//
// The simplest interpolator, and the place to meet the two conventions the
// whole library shares: construction takes ranges, and evaluation is
// Evaluate<N> with operator() as a shorthand for Evaluate<0>.

#include <Interpolation/Linear.hpp>

#include <cstdio>
#include <stdexcept>
#include <vector>

int
main() {
    // Abscissae must be strictly increasing. They need not be evenly spaced,
    // and here they deliberately are not.
    const std::vector<double> x{0.0, 0.5, 2.0, 3.5, 4.0};
    const std::vector<double> y{0.0, 1.0, 3.0, 2.0, 0.0};

    const Interpolation::Linear f{x, y};

    std::printf("A piecewise-linear interpolant through %zu points\n\n",
                f.Size());
    std::printf("      x        f(x)      f'(x)\n");
    for (int i = 0; i <= 8; ++i) {
        const double q = 4.0 * i / 8.0;
        std::printf("  %7.4f  %9.5f  %9.5f\n", q, f(q), f.Evaluate<1>(q));
    }

    // The interpolant is linear on each segment, so the second derivative and
    // everything above it are identically zero.
    std::printf("\n  f''(1.0) = %.1f, and every higher order is zero too\n",
                f.Evaluate<2>(1.0));

    // Queries outside the sample range continue the first or last segment
    // rather than failing. That is extrapolation, and it is on you to decide
    // whether it means anything for your data.
    std::printf("\n  Outside the samples the end segment continues:\n");
    std::printf("    f(-1.0) = %8.5f     f(5.0) = %8.5f\n", f(-1.0), f(5.0));

    // Input is checked when the object is built, not when it is used, and a
    // problem is reported by throwing rather than by an assert that would
    // vanish in a release build.
    std::printf("\n  Bad input is rejected at construction:\n");
    try {
        const std::vector<double> repeated{0.0, 1.0, 1.0, 2.0};
        const std::vector<double> values{0.0, 1.0, 2.0, 3.0};
        const Interpolation::Linear bad{repeated, values};
        (void) bad;
    } catch (const std::invalid_argument &error) {
        std::printf("    %s\n", error.what());
    }

    return 0;
}
