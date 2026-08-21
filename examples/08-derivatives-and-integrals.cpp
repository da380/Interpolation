// 08 - Derivatives and integrals
//
// Derivative and Primitive both turn a function into another function, so an
// antiderivative can be differentiated back, added to something else, or
// integrated over a range.

#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Function.hpp>
#include <Interpolation/Linear.hpp>
#include <Interpolation/Polynomial.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

int
main() {
    using namespace Interpolation;

    // Sample a cubic, and clamp the spline with the true end slopes so that
    // it reproduces the cubic exactly. Then its integral is exact too.
    const auto cubic = [](double v) { return v * v * v; };
    const auto slope = [](double v) { return 3.0 * v * v; };
    const std::vector<double> x{0.0, 0.7, 1.9, 2.4, 3.0};
    std::vector<double> y;
    for (const auto v : x) {
        y.push_back(cubic(v));
    }

    const CubicSpline s{x, y, BoundaryCondition::Clamped, slope(0.0),
                        slope(3.0)};

    // Primitive is the antiderivative vanishing at the first node. The
    // cumulative integral is worked out once, when the node is built, so an
    // interpolator carries none of that unless it is asked for.
    const auto area = Primitive(s);

    std::printf("Integrating a spline through data sampled from x^3\n\n");
    std::printf("  integral over [0, 3] = %.10f\n", area.Integral(0.0, 3.0));
    std::printf("  exact (81/4)         = %.10f\n", 81.0 / 4.0);

    // The antiderivative is itself a function, so differentiating it returns
    // the original.
    const auto back = Derivative(area);
    std::printf("\n  d/dx of the antiderivative, against the spline:\n");
    for (const double q : {0.5, 1.5, 2.5}) {
        std::printf("    x = %.1f   %.10f   %.10f\n", q, back(q), s(q));
    }

    // The piecewise-linear interpolant integrates to the trapezoid rule, by
    // construction rather than by coincidence.
    const Linear l{x, y};
    const auto trapezoidal = Primitive(l);
    double byHand = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        byHand += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2.0;
    }
    std::printf("\n  Linear interpolant over [0, 3]: %.10f\n",
                trapezoidal.Integral(0.0, 3.0));
    std::printf("  trapezoid rule by hand:         %.10f\n", byHand);

    // A polynomial has a closed-form antiderivative, so it takes a different
    // route to the same interface.
    const Polynomial<double> p{0.0, 0.0, 0.0, 1.0}; // x^3
    std::printf("\n  Polynomial x^3 over [0, 3]:     %.10f\n",
                Primitive(p).Integral(0.0, 3.0));

    return 0;
}
