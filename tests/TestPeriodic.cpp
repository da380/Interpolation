#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "AllocationCounter.hpp"
#include "TestCubicSpline.h"
#include "TestUtilities.h"

// The periodic end condition. The final node closes the period, so a grid
// the caller already holds needs no special case; the system is cyclic and is
// solved by Sherman-Morrison on the one tridiagonal factorisation; queries
// are reduced into one period, and everything that answers a query — the
// segment search, the nodal sweep, the antiderivative — has to agree about
// what that means.

namespace {

using Interpolation::BoundaryCondition;
using Interpolation::CubicSpline;
using Interpolation::CubicSplineSystem;
using Interpolation::Side;

constexpr double kPeriod = 2.0 * std::numbers::pi;

// One period on deliberately uneven nodes, the closing node included.
const std::vector<double> kNodes{0.0, 0.5, 1.1, 1.3,    2.4,
                                 3.9, 4.2, 5.0, kPeriod};

// Smooth and periodic, with more than one harmonic so that no curvature is
// accidentally zero.
double
Wave(double x) {
    return std::sin(x) + 0.5 * std::cos(2.0 * x) + 0.2 * std::sin(3.0 * x);
}

// Samples of `f` on `x`, with the closing ordinate copied from the first
// rather than evaluated, which is the contract.
template <typename F>
auto
Closed(const std::vector<double> &x, F f) {
    std::vector<decltype(f(0.0))> y;
    y.reserve(x.size());
    for (const auto value : x) {
        y.push_back(f(value));
    }
    y.back() = y.front();
    return y;
}

// `n` segments over one period: uniform, or stretched by the monotone map
// t + 0.4 sin t, which keeps the spacing uneven at every resolution without
// making it irregular, so a convergence rate still means something.
std::vector<double>
PeriodGrid(std::size_t n, bool stretched) {
    std::vector<double> x(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        const auto t =
            kPeriod * static_cast<double>(i) / static_cast<double>(n);
        x[i] = stretched ? t + 0.4 * std::sin(t) : t;
    }
    x.front() = 0.0;
    x.back() = kPeriod;
    return x;
}

CubicSpline<std::ranges::ref_view<const std::vector<double>>,
            std::ranges::ref_view<const std::vector<double>>>
PeriodicSpline(const std::vector<double> &x, const std::vector<double> &y) {
    return CubicSpline{x, y, BoundaryCondition::Periodic, 0.0, 0.0};
}

} // namespace

TEST(PeriodicSpline, CurvaturesAgreeWithADenseSolveOfTheCyclicSystem) {
    const auto y = Closed(kNodes, Wave);
    const CubicSplineSystem system{kNodes, BoundaryCondition::Periodic};
    const auto curvature = system.Curvatures(y);

    // The cyclic system written out in full, corners and all, and solved by
    // pivoted elimination: nothing here knows it is nearly tridiagonal.
    const auto m = kNodes.size() - 1;
    std::vector<double> h(m), secant(m);
    for (std::size_t i = 0; i < m; ++i) {
        h[i] = kNodes[i + 1] - kNodes[i];
        secant[i] = (y[i + 1] - y[i]) / h[i];
    }
    std::vector<std::vector<double>> a(m, std::vector<double>(m, 0.0));
    std::vector<double> b(m);
    for (std::size_t i = 0; i < m; ++i) {
        const auto previous = (i + m - 1) % m;
        a[i][previous] += h[previous] / 6.0;
        a[i][i] += (h[previous] + h[i]) / 3.0;
        a[i][(i + 1) % m] += h[i] / 6.0;
        b[i] = secant[i] - secant[previous];
    }
    const auto reference = CubicSplineTest::SolveDense(a, b);

    ASSERT_EQ(curvature.size(), kNodes.size());
    for (std::size_t i = 0; i < m; ++i) {
        SCOPED_TRACE(i);
        InterpolationTest::ExpectScaledNear(curvature[i], reference[i]);
    }
    EXPECT_EQ(curvature[m], curvature[0])
        << "the closing node's curvature is the first's, by construction";
}

TEST(PeriodicSpline, ConvergesAtFourthOrder) {
    // Measured worst errors on a fine sample, halving the spacing each time.
    // Rates observed for sin and cos, on both grids: 4.00 to 4.10.
    for (const bool stretched : {false, true}) {
        for (const auto f : {+[](double x) { return std::sin(x); },
                             +[](double x) { return std::cos(x); }}) {
            std::vector<double> worst;
            for (const std::size_t n : {16, 32, 64, 128}) {
                const auto x = PeriodGrid(n, stretched);
                const auto y = Closed(x, f);
                const auto spline = PeriodicSpline(x, y);
                double error = 0.0;
                for (int k = 0; k <= 4000; ++k) {
                    const auto q = kPeriod * k / 4000.0;
                    error = std::max(error, std::abs(spline(q) - f(q)));
                }
                worst.push_back(error);
            }
            for (std::size_t i = 0; i + 1 < worst.size(); ++i) {
                SCOPED_TRACE(std::string{stretched ? "stretched" : "uniform"} +
                             ", refinement " + std::to_string(i));
                const auto rate = std::log2(worst[i] / worst[i + 1]);
                EXPECT_GT(rate, 3.8);
                EXPECT_LT(rate, 4.2);
            }
        }
    }
}

TEST(PeriodicSpline, IsTwiceContinuouslyDifferentiableAcrossTheSeam) {
    const auto y = Closed(kNodes, Wave);
    const auto spline = PeriodicSpline(kNodes, y);

    // From either side of the seam, through the queries themselves: the last
    // piece read just below the closing node against the first piece read
    // just above the first node. A not-a-knot spline through the same data
    // misses by 0.17 in slope and 6.4 in curvature.
    const auto delta = 1.0e-10;
    const auto below = kPeriod - delta;
    EXPECT_NEAR(spline.Evaluate<0>(below), spline.Evaluate<0>(delta), 1.0e-8);
    EXPECT_NEAR(spline.Evaluate<1>(below), spline.Evaluate<1>(delta), 1.0e-8);
    EXPECT_NEAR(spline.Evaluate<2>(below), spline.Evaluate<2>(delta), 1.0e-8);

    // And exactly at the seam, through the nodal sweep, where the left limit
    // at the first node is the last piece's right end.
    for (const auto side : {Side::Left, Side::Right}) {
        SCOPED_TRACE(side == Side::Left ? "left" : "right");
        const auto values = spline.NodeValues<0>(side);
        const auto slopes = spline.NodeValues<1>(side);
        const auto curvatures = spline.NodeValues<2>(side);
        InterpolationTest::ExpectScaledNear(values.front(), values.back());
        InterpolationTest::ExpectScaledNear(slopes.front(), slopes.back());
        InterpolationTest::ExpectScaledNear(curvatures.front(),
                                            curvatures.back());
    }
    const auto left = spline.NodeValues<1>(Side::Left);
    const auto right = spline.NodeValues<1>(Side::Right);
    InterpolationTest::ExpectScaledNear(left.front(), right.front());
}

TEST(PeriodicSpline, AQueryMovedByWholePeriodsGivesTheSameAnswer) {
    const auto y = Closed(kNodes, Wave);
    const auto spline = PeriodicSpline(kNodes, y);

    // Queries kept clear of the knots, so the third derivative, which is
    // piecewise constant, is not asked to choose between two pieces.
    for (const double q : {0.2, 1.2, 3.0, 4.6, 6.1}) {
        for (const int turns : {-3, -1, 1, 2, 5}) {
            SCOPED_TRACE(std::to_string(q) + " moved by " +
                         std::to_string(turns));
            const auto moved = q + turns * kPeriod;
            EXPECT_EQ(spline.Segment(moved), spline.Segment(q));
            InterpolationTest::ExpectScaledNear(spline.Evaluate<0>(moved),
                                                spline.Evaluate<0>(q));
            InterpolationTest::ExpectScaledNear(spline.Evaluate<1>(moved),
                                                spline.Evaluate<1>(q));
            InterpolationTest::ExpectScaledNear(spline.Evaluate<2>(moved),
                                                spline.Evaluate<2>(q));
            InterpolationTest::ExpectScaledNear(spline.Evaluate<3>(moved),
                                                spline.Evaluate<3>(q));
        }
    }

    // The period is half-open: the closing node is the first node again.
    EXPECT_EQ(spline.Segment(kPeriod), 0u);
    EXPECT_EQ(spline(kPeriod), y.front());
}

TEST(PeriodicSpline, TheNodalSweepAgreesWithTheQueryAtEveryNode) {
    const auto y = Closed(kNodes, Wave);
    const auto spline = PeriodicSpline(kNodes, y);
    const auto n = kNodes.size();

    EXPECT_EQ(spline.NodeValues<0>(), y);

    // Right is what a query at a node selects, and at the closing node a
    // query is reduced to the first, so the sweep must wrap there too.
    const auto third = spline.NodeValues<3>(Side::Right);
    for (std::size_t k = 0; k < n; ++k) {
        SCOPED_TRACE(k);
        EXPECT_EQ(third[k], spline.Evaluate<3>(spline.Node(k)));
    }

    // The two limits of the third derivative at the seam are those of the
    // first and the last piece, read anywhere inside them.
    const auto firstPiece = spline.Evaluate<3>(0.5 * kNodes[1]);
    const auto lastPiece =
        spline.Evaluate<3>(0.5 * (kNodes[n - 2] + kNodes[n - 1]));
    ASSERT_NE(firstPiece, lastPiece) << "the seam is a knot like any other";
    EXPECT_EQ(third[n - 1], firstPiece);
    EXPECT_EQ(spline.NodeValues<3>(Side::Left)[0], lastPiece);
}

TEST(PeriodicSpline, WorksForComplexOrdinates) {
    const auto real = Closed(kNodes, Wave);
    const auto imaginary =
        Closed(kNodes, [](double x) { return std::cos(x) - std::sin(2 * x); });
    std::vector<std::complex<double>> y(kNodes.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        y[i] = {real[i], imaginary[i]};
    }

    const CubicSpline spline{kNodes, y, BoundaryCondition::Periodic,
                             std::complex<double>{}, std::complex<double>{}};
    const auto realSpline = PeriodicSpline(kNodes, real);
    const auto imaginarySpline = PeriodicSpline(kNodes, imaginary);

    // The matrix is real, so a complex spline is two real ones side by side.
    for (const double q : {-1.0, 0.3, 2.2, 4.4, 5.9, 9.0}) {
        SCOPED_TRACE(q);
        InterpolationTest::ExpectScaledNear(
            spline.Evaluate<1>(q),
            std::complex<double>{realSpline.Evaluate<1>(q),
                                 imaginarySpline.Evaluate<1>(q)});
    }
}

TEST(PeriodicSpline, SolvingAllocatesNothing) {
    if constexpr (!InterpolationTest::AllocationCountingEnabled()) {
        GTEST_SKIP() << "allocation counting stands down under sanitizers";
    } else {
        const CubicSplineSystem system{kNodes, BoundaryCondition::Periodic};
        const auto y = Closed(kNodes, Wave);
        std::vector<double> curvature(kNodes.size());
        std::vector<double> nodal(kNodes.size());

        // The correction is solved for once, at construction, and kept.
        const auto before = InterpolationTest::AllocationCount();
        for (int line = 0; line < 64; ++line) {
            system.Solve(y, std::span<double>{curvature});
            system.EvaluateAtNodes<3>(y, curvature, std::span<double>{nodal},
                                      Side::Left);
        }
        const auto after = InterpolationTest::AllocationCount();

        EXPECT_EQ(after, before)
            << "solving allocated " << (after - before) << " times";
    }
}

TEST(PeriodicSpline, ThePrimitiveAdvancesByOnePeriodIntegralPerPeriod) {
    using Interpolation::Derivative;
    using Interpolation::Primitive;
    using Interpolation::Ref;

    // A mean of one, so the primitive grows by about 2 pi a period and an
    // answer that wrapped instead of growing could not pass.
    const auto offset = [](double x) { return 1.0 + Wave(x); };
    const auto y = Closed(kNodes, offset);
    const auto spline = PeriodicSpline(kNodes, y);
    const auto S = Primitive(spline);

    const auto onePeriod = S(kPeriod);
    EXPECT_NEAR(onePeriod, kPeriod, 0.05);

    for (const double q : {0.2, 1.2, 3.0, 4.6, 6.1}) {
        for (const int turns : {-2, -1, 1, 3}) {
            SCOPED_TRACE(std::to_string(q) + " moved by " +
                         std::to_string(turns));
            const auto moved = q + turns * kPeriod;
            EXPECT_NEAR(S(moved) - S(q), turns * onePeriod, 1.0e-12);
            // It still differentiates back to the spline out there.
            InterpolationTest::ExpectScaledNear(Derivative(S)(moved),
                                                spline(q));
            // Any window one period wide holds the same integral.
            EXPECT_NEAR(S.Integral(moved, moved + kPeriod), onePeriod, 1.0e-12);
        }
    }

    // Borrowing the spline must not lose the period on the way in.
    const auto borrowed = Primitive(Ref(spline));
    for (const double q : {-7.0, 2.5, 20.0}) {
        SCOPED_TRACE(q);
        EXPECT_EQ(borrowed(q), S(q));
    }
}

TEST(PeriodicSpline, ReportsItsPeriod) {
    const CubicSplineSystem periodic{kNodes, BoundaryCondition::Periodic};
    ASSERT_TRUE(periodic.Period().has_value());
    EXPECT_EQ(*periodic.Period(), kNodes.back() - kNodes.front());

    const CubicSplineSystem natural{kNodes, BoundaryCondition::Natural};
    EXPECT_FALSE(natural.Period().has_value());

    const auto y = Closed(kNodes, Wave);
    EXPECT_EQ(PeriodicSpline(kNodes, y).Period(), periodic.Period());
}

// Refused, and in the words of the class that was called.
template <typename Make>
void
ExpectRefusedBy(const char *name, Make make) {
    try {
        make();
        ADD_FAILURE() << "expected std::invalid_argument";
    } catch (const std::invalid_argument &error) {
        EXPECT_TRUE(std::string_view{error.what()}.starts_with(name))
            << error.what();
    }
}

TEST(PeriodicSpline, RefusesWhatItCannotUse) {
    const auto y = Closed(kNodes, Wave);

    // A closing ordinate that is close is not a copy.
    auto open = y;
    open.back() = std::nextafter(open.back(), 10.0);
    ExpectRefusedBy("CubicSpline:",
                    [&] { (void) PeriodicSpline(kNodes, open); });
    const CubicSplineSystem system{kNodes, BoundaryCondition::Periodic};
    std::vector<double> curvature(kNodes.size());
    ExpectRefusedBy("CubicSplineSystem::Solve:",
                    [&] { system.Solve(open, std::span<double>{curvature}); });

    // Three nodes, the closing one included, is two distinct points.
    const std::vector<double> few{0.0, 1.0, 2.0};
    const std::vector<double> fewValues{1.0, 2.0, 1.0};
    ExpectRefusedBy("CubicSpline:",
                    [&] { (void) PeriodicSpline(few, fewValues); });
    ExpectRefusedBy("CubicSplineSystem:", [&] {
        (void) CubicSplineSystem{few, BoundaryCondition::Periodic};
    });

    // Periodic is a condition on the whole system, so not at one end only.
    for (const auto other :
         {BoundaryCondition::Natural, BoundaryCondition::Clamped,
          BoundaryCondition::NotAKnot}) {
        ExpectRefusedBy("CubicSpline:", [&] {
            (void) CubicSpline{kNodes, y,     BoundaryCondition::Periodic,
                               0.0,    other, 0.0};
        });
        ExpectRefusedBy("CubicSpline:", [&] {
            (void) CubicSpline{
                kNodes, y, other, 0.0, BoundaryCondition::Periodic, 0.0};
        });
        ExpectRefusedBy("CubicSplineSystem:", [&] {
            (void) CubicSplineSystem{kNodes, BoundaryCondition::Periodic,
                                     other};
        });
    }
}
