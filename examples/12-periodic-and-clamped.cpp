// 12 - Periodic and clamped
//
// Data on a circle has no ends, and a spline through it should not invent
// them. BoundaryCondition::Periodic joins the final node to the first with
// the value and two derivatives continuous across the seam. On a grid, each
// axis takes its own conditions: here a field over colatitude and longitude,
// clamped at the poles with a different slope on every meridian, periodic in
// longitude.

#include <Interpolation/BicubicSpline.hpp>
#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Function.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <utility>
#include <vector>

int
main() {
    using namespace Interpolation;
    constexpr double pi = std::numbers::pi;

    // The final node closes the period, one period after the first, and its
    // ordinate is a copy of the first: copied, not recomputed, since a
    // periodic spline refuses a closing value that differs by rounding.
    const auto wave = [](double x) {
        return 1.0 + std::sin(x) + 0.5 * std::cos(2 * x);
    };
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i <= 12; ++i) {
        x.push_back(2.0 * pi * i / 12.0);
        y.push_back(wave(x.back()));
    }
    y.back() = y.front();

    const CubicSpline periodic{
        x,  y, BoundaryCondition::Periodic, 0.0, BoundaryCondition::Periodic,
        0.0};
    const CubicSpline notAKnot{
        x,  y, BoundaryCondition::NotAKnot, 0.0, BoundaryCondition::NotAKnot,
        0.0};

    // The slope arriving at 2 pi along the last piece, and leaving 0 along
    // the first. For the periodic spline those are the two sides of one
    // point; not-a-knot never heard of the seam.
    std::printf("Slope arriving at 2 pi, and leaving 0:\n");
    for (const auto &[name, spline] : {std::pair{"periodic", &periodic},
                                       std::pair{"not-a-knot", &notAKnot}}) {
        std::printf("  %-10s  %9.6f  %9.6f\n", name,
                    spline->NodeValues<1>(Side::Left).back(),
                    spline->NodeValues<1>(Side::Right).front());
    }

    // Queries are reduced into one period, so the spline repeats ...
    std::printf("\nThe same point, whole periods apart:\n");
    for (const double q : {1.0, 1.0 + 2.0 * pi, 1.0 - 6.0 * pi}) {
        std::printf("  s(%8.4f) = %.12f\n", q, periodic(q));
    }

    // ... and its antiderivative does not: it gains one period's integral per
    // period, as the antiderivative of a function with nonzero mean must.
    const auto integral = Primitive(periodic);
    std::printf("\nIntegral over one period: %.8f (exact %.8f)\n",
                integral(2.0 * pi), 2.0 * pi);
    std::printf("Integral over [1, 1 + 3 period]: %.8f\n",
                integral.Integral(1.0, 1.0 + 6.0 * pi));

    // A field over colatitude theta in [0, pi] and longitude phi in [0, 2 pi].
    const auto field = [](double theta, double phi) {
        return std::cos(theta) + 0.3 * std::sin(theta) * std::cos(phi);
    };
    const auto fieldDTheta = [](double theta, double phi) {
        return -std::sin(theta) + 0.3 * std::cos(theta) * std::cos(phi);
    };

    std::vector<double> theta;
    std::vector<double> phi;
    for (int i = 0; i <= 16; ++i) {
        theta.push_back(pi * i / 16.0);
    }
    for (int j = 0; j <= 32; ++j) {
        phi.push_back(2.0 * pi * j / 32.0);
    }

    // Values row-major, the final column a copy of the first. The pole slopes
    // are d/dtheta at each meridian, and they close in longitude just as the
    // values do.
    std::vector<double> values;
    for (const auto t : theta) {
        for (const auto p : phi) {
            values.push_back(field(t, p));
        }
        values.back() = values[values.size() - phi.size()];
    }
    std::vector<double> north;
    std::vector<double> south;
    for (const auto p : phi) {
        north.push_back(fieldDTheta(0.0, p));
        south.push_back(fieldDTheta(pi, p));
    }
    north.back() = north.front();
    south.back() = south.front();

    const BicubicSpline sphere{
        theta,
        phi,
        values,
        AxisConditions{BoundaryCondition::Clamped, BoundaryCondition::Clamped},
        AxisConditions{BoundaryCondition::Periodic,
                       BoundaryCondition::Periodic},
        EdgeSlopes<double>{.xLow = north, .xHigh = south}};

    std::printf("\nOn a %zu x %zu (theta, phi) grid, against exact:\n",
                theta.size(), phi.size());
    for (const auto &[t, p] :
         {std::pair{0.7, 0.01}, std::pair{0.7, 2.0 * pi - 0.01},
          std::pair{0.7, 0.01 + 2.0 * pi}, std::pair{0.05, 3.0}}) {
        std::printf("  (%.2f, %7.4f)  %10.7f  error %.1e\n", t, p, sphere(t, p),
                    sphere(t, p) - field(t, p));
    }
    std::printf("  d/dtheta at the north pole, between meridians: %.7f "
                "(exact %.7f)\n",
                sphere.Evaluate<1, 0>(0.0, 0.1), fieldDTheta(0.0, 0.1));
}
