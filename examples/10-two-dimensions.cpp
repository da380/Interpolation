// 10 - Two dimensions
//
// Bilinear and BicubicSpline interpolate on a rectilinear grid: two axes and
// a flat, row-major value range. Both are tensor products of the
// corresponding one-dimensional scheme, so a mixed partial derivative is just
// the two one-dimensional rules applied in turn.

#include <Interpolation/BicubicSpline.hpp>
#include <Interpolation/Bilinear.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int
main() {
    using namespace Interpolation;

    const auto exact = [](double a, double b) {
        return std::sin(a) * std::cos(b);
    };

    // Axes need not be evenly spaced, and the two need not match.
    const int nx = 9;
    const int ny = 9;
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < nx; ++i) {
        x.push_back(3.0 * i / (nx - 1));
    }
    for (int j = 0; j < ny; ++j) {
        y.push_back(3.0 * j / (ny - 1));
    }

    // Values are flat and row-major: element (i, j) sits at i * size(y) + j.
    std::vector<double> v;
    v.reserve(x.size() * y.size());
    for (const auto a : x) {
        for (const auto b : y) {
            v.push_back(exact(a, b));
        }
    }

    const Bilinear bilinear{x, y, v};
    const BicubicSpline bicubic{x, y, v};

    std::printf("sin(x)cos(y) on a %d x %d grid\n\n", nx, ny);
    std::printf("     x      y       exact    bilinear       error     bicubic"
                "       error\n");
    double worstBilinear = 0.0;
    double worstBicubic = 0.0;
    for (int i = 0; i <= 5; ++i) {
        const double px = 0.4 + 2.2 * i / 5.0;
        const double py = 0.7 + 1.5 * i / 5.0;
        const double truth = exact(px, py);
        worstBilinear =
            std::max(worstBilinear, std::abs(bilinear(px, py) - truth));
        worstBicubic =
            std::max(worstBicubic, std::abs(bicubic(px, py) - truth));
        std::printf("  %5.3f  %5.3f  %10.6f  %10.6f  %10.2e  %10.6f  %10.2e\n",
                    px, py, truth, bilinear(px, py), bilinear(px, py) - truth,
                    bicubic(px, py), bicubic(px, py) - truth);
    }
    std::printf("\n  worst here: bilinear %.2e, bicubic %.2e\n", worstBilinear,
                worstBicubic);

    // Mixed partial derivatives come from the same object.
    const double px = 1.2;
    const double py = 0.8;
    std::printf(
        "\nDerivatives of the bicubic at (%.1f, %.1f), against exact:\n", px,
        py);
    std::printf("  d/dx    %10.6f  (%10.6f)\n", bicubic.Evaluate<1, 0>(px, py),
                std::cos(px) * std::cos(py));
    std::printf("  d/dy    %10.6f  (%10.6f)\n", bicubic.Evaluate<0, 1>(px, py),
                -std::sin(px) * std::sin(py));
    std::printf("  d2/dxdy %10.6f  (%10.6f)\n", bicubic.Evaluate<1, 1>(px, py),
                -std::cos(px) * std::sin(py));

    // The edge condition matters. Not-a-knot is the default because it keeps
    // the scheme fourth order at the boundary, where the natural condition
    // forces a curvature the data does not have and costs an order.
    const BicubicSpline natural{x, y, v, BoundaryCondition::Natural};
    double edgeNatural = 0.0;
    double edgeNotAKnot = 0.0;
    for (int i = 0; i <= 40; ++i) {
        const double q = 3.0 * i / 40.0;
        const double truth = exact(q, 0.0);
        edgeNatural = std::max(edgeNatural, std::abs(natural(q, 0.0) - truth));
        edgeNotAKnot =
            std::max(edgeNotAKnot, std::abs(bicubic(q, 0.0) - truth));
    }
    std::printf("\nWorst error along the y = 0 edge:\n");
    std::printf("  natural    %.3e\n  not-a-knot %.3e\n", edgeNatural,
                edgeNotAKnot);

    return 0;
}
