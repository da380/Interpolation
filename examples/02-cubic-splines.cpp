// 02 - Cubic splines
//
// A cubic spline is smooth where the linear interpolant has corners: value,
// slope and curvature are all continuous across the nodes. That is what makes
// it the usual choice for a physical profile.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Linear.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int
main() {
    // Sample a smooth function on a coarse, uneven grid.
    const auto exact = [](double v) { return std::sin(v); };
    const std::vector<double> x{0.0, 0.9, 1.7, 2.2, 3.4, 4.1, 5.0, 6.0};
    std::vector<double> y;
    for (const auto v : x) {
        y.push_back(exact(v));
    }

    const Interpolation::CubicSpline spline{x, y};
    const Interpolation::Linear linear{x, y};

    std::printf("Interpolating sin(x) from %zu uneven samples\n\n",
                spline.Size());
    std::printf("      x       exact      spline       error      linear "
                "      error\n");
    double worstSpline = 0.0;
    double worstLinear = 0.0;
    for (int i = 0; i <= 10; ++i) {
        const double q = 6.0 * i / 10.0;
        const double truth = exact(q);
        const double s = spline(q);
        const double l = linear(q);
        worstSpline = std::max(worstSpline, std::abs(s - truth));
        worstLinear = std::max(worstLinear, std::abs(l - truth));
        std::printf("  %6.3f  %10.6f  %10.6f  %10.2e  %10.6f  %10.2e\n", q,
                    truth, s, s - truth, l, l - truth);
    }
    std::printf("\n  worst error: spline %.2e, linear %.2e\n", worstSpline,
                worstLinear);

    // Evaluate<N> gives derivatives from the same object. The pieces are
    // cubic, so the third derivative is piecewise constant and the fourth
    // vanishes.
    std::printf("\n  At x = 2.0, against the exact derivatives of sin:\n");
    std::printf("    f    %10.6f   (%10.6f)\n", spline(2.0), std::sin(2.0));
    std::printf("    f'   %10.6f   (%10.6f)\n", spline.Evaluate<1>(2.0),
                std::cos(2.0));
    std::printf("    f''  %10.6f   (%10.6f)\n", spline.Evaluate<2>(2.0),
                -std::sin(2.0));
    std::printf("    f''' %10.6f   (piecewise constant)\n",
                spline.Evaluate<3>(2.0));
    std::printf("    f''''%10.6f   (identically zero)\n",
                spline.Evaluate<4>(2.0));

    return 0;
}
