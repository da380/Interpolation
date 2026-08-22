// 04 - Akima splines and Lagrange interpolation
//
// Two interpolants with quite different characters: one local and resistant
// to overshoot, one global and exact for polynomials.

#include <Interpolation/AkimaSpline.hpp>
#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Lagrange.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int
main() {
    // Data with a flat shelf and a step. A cubic spline rings on this; Akima
    // weights its slopes by local secant differences and largely does not.
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const std::vector<double> y{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    const Interpolation::AkimaSpline akima{x, y};
    const Interpolation::CubicSpline spline{x, y};

    std::printf("A step in the data: Akima against a cubic spline\n\n");
    std::printf("      x       akima      spline\n");
    double akimaUnder = 0.0;
    double splineUnder = 0.0;
    for (int i = 0; i <= 14; ++i) {
        const double q = 7.0 * i / 14.0;
        akimaUnder = std::min(akimaUnder, akima(q));
        splineUnder = std::min(splineUnder, spline(q));
        std::printf("  %6.3f  %10.6f  %10.6f\n", q, akima(q), spline(q));
    }
    std::printf("\n  most negative value reached: akima %.6f, spline %.6f\n",
                akimaUnder, splineUnder);
    std::printf("  The data never goes below zero; overshoot below it is an\n"
                "  artefact of the interpolant, not a feature of the data.\n");

    // Lagrange builds the single polynomial through every node, so with n
    // nodes it reproduces any polynomial of degree n-1 exactly.
    const std::vector<double> nodes{-1.0, 0.25, 2.0, 4.0};
    const auto cubic = [](double v) {
        return 1.0 - 2.0 * v + 0.5 * v * v + 1.25 * v * v * v;
    };
    std::vector<double> values;
    for (const auto v : nodes) {
        values.push_back(cubic(v));
    }
    const Interpolation::Lagrange global{nodes, values};

    std::printf("\nLagrange through 4 nodes, on data from a cubic\n\n");
    std::printf("      x        exact  interpolant       error\n");
    for (const double q : {-2.0, -0.5, 1.0, 3.0, 5.0}) {
        std::printf("  %6.3f  %11.5f  %11.5f  %10.2e\n", q, cubic(q), global(q),
                    global(q) - cubic(q));
    }
    std::printf("\n  Exact everywhere, including outside the nodes, because a\n"
                "  cubic through four points is the cubic. That stops being\n"
                "  comfortable at high degree, where the polynomial through\n"
                "  many points oscillates badly between them.\n");

    // The cardinal basis is available on its own. Basis function i is one at
    // node i and zero at every other node, and the set sums to one.
    const Interpolation::LagrangeBasis basis{nodes};
    std::printf("\nCardinal basis at x = 1.0:\n   ");
    double total = 0.0;
    for (std::size_t i = 0; i < basis.Size(); ++i) {
        total += basis(i, 1.0);
        std::printf(" L%zu = %8.5f", i, basis(i, 1.0));
    }
    std::printf("\n    sum = %.15f\n", total);

    return 0;
}
