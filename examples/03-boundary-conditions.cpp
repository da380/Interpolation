// 03 - Boundary conditions
//
// A cubic spline needs two more conditions than the data supplies, and the
// choice matters near the ends. This compares the three on offer against data
// that is exactly cubic, where the right answer is known.

#include <Interpolation/CubicSpline.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int
main() {
    using Interpolation::BoundaryCondition;
    using Interpolation::CubicSpline;

    // A genuine cubic, so a well-chosen spline should reproduce it exactly.
    const auto exact = [](double v) {
        return 2.0 - 0.5 * v + 0.25 * v * v + 1.5 * v * v * v;
    };
    const auto slope = [](double v) { return -0.5 + 0.5 * v + 4.5 * v * v; };

    const std::vector<double> x{0.0, 0.4, 1.3, 1.9, 3.0, 4.2, 5.0};
    std::vector<double> y;
    for (const auto v : x) {
        y.push_back(exact(v));
    }

    // Natural: the second derivative is forced to zero at both ends. Cheap,
    // and wrong here, because the data has curvature there.
    const CubicSpline natural{x, y};

    // Clamped: you supply the end slopes. Exact when you know them.
    const CubicSpline clamped{x, y, BoundaryCondition::Clamped, slope(0.0),
                              slope(5.0)};

    // Not-a-knot: the third derivative is made continuous across the first and
    // last interior knots, so nothing false is imposed. It needs no extra
    // information and needs at least four nodes.
    const CubicSpline notAKnot{
        x,  y, BoundaryCondition::NotAKnot, 0.0, BoundaryCondition::NotAKnot,
        0.0};

    std::printf(
        "Data sampled from an exact cubic; worst error over [0, 5]\n\n");
    std::printf("      x        exact     natural     clamped  not-a-knot\n");
    double wn = 0.0;
    double wc = 0.0;
    double wk = 0.0;
    for (int i = 0; i <= 10; ++i) {
        const double q = 5.0 * i / 10.0;
        const double truth = exact(q);
        wn = std::max(wn, std::abs(natural(q) - truth));
        wc = std::max(wc, std::abs(clamped(q) - truth));
        wk = std::max(wk, std::abs(notAKnot(q) - truth));
        std::printf("  %6.3f  %11.5f  %10.2e  %10.2e  %10.2e\n", q, truth,
                    natural(q) - truth, clamped(q) - truth,
                    notAKnot(q) - truth);
    }
    std::printf(
        "\n  worst:  natural %8.2e   clamped %8.2e   not-a-knot %8.2e\n", wn,
        wc, wk);
    std::printf(
        "\n  Clamped and not-a-knot recover the cubic exactly. Natural\n"
        "  cannot: it insists the curvature vanishes at the ends, which\n"
        "  for this data is simply untrue, and the error spreads inward.\n");

    return 0;
}
