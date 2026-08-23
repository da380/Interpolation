#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <cmath>
#include <complex>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "AllocationCounter.hpp"
#include "TestUtilities.h"

// Evaluating at every node is what a differentiation operator on a fixed grid
// actually asks for, and it is the one query that needs no search: the segment
// adjoining each node is known. These checks pin the path to the general one
// it shortcuts, so the two cannot drift apart.

namespace {

using Interpolation::AkimaSpline;
using Interpolation::BoundaryCondition;
using Interpolation::CubicSpline;
using Interpolation::CubicSplineSystem;
using Interpolation::Linear;
using Interpolation::Side;

// Uneven spacing throughout, so nothing here can pass by accident on a
// uniform grid.
const std::vector<double> kNodes{0.0, 0.4, 1.5, 1.7, 2.9, 4.4, 4.6, 6.0};

std::vector<double>
SampledOn(const std::vector<double> &x) {
    std::vector<double> y;
    y.reserve(x.size());
    for (const auto value : x) {
        y.push_back(std::exp(-0.3 * value) * std::sin(1.7 * value));
    }
    return y;
}

// Somewhere strictly inside segment i, where a piecewise-constant derivative
// takes its value for that segment.
double
Inside(std::size_t i) {
    return 0.5 * (kNodes[i] + kNodes[i + 1]);
}

// The nodal sweep must agree with the general query at every node, for every
// derivative order the interpolant carries.
template <std::size_t N, typename Interpolant>
void
CheckOrderAgrees(const Interpolant &f) {
    SCOPED_TRACE("derivative order " + std::to_string(N));
    const auto n = f.Size();

    std::vector<typename Interpolant::Scalar> swept(n);
    f.template EvaluateAtNodes<N>(
        std::span<typename Interpolant::Scalar>{swept});

    for (std::size_t k = 0; k < n; ++k) {
        SCOPED_TRACE("node " + std::to_string(k));
        // Bit-for-bit. The two paths run the same segment formula on the same
        // segment; anything less would mean one of them had grown its own.
        EXPECT_EQ(swept[k], f.template Evaluate<N>(f.Node(k)));
    }

    const auto returned = f.template NodeValues<N>();
    ASSERT_EQ(returned.size(), n);
    for (std::size_t k = 0; k < n; ++k) {
        EXPECT_EQ(returned[k], swept[k]);
    }
}

} // namespace

TEST(EvaluateAtNodes, AgreesWithTheGeneralQueryForEveryInterpolant) {
    const auto y = SampledOn(kNodes);

    {
        const Linear f{kNodes, y};
        CheckOrderAgrees<0>(f);
        CheckOrderAgrees<1>(f);
        CheckOrderAgrees<2>(f);
    }
    {
        const CubicSpline f{kNodes, y};
        CheckOrderAgrees<0>(f);
        CheckOrderAgrees<1>(f);
        CheckOrderAgrees<2>(f);
        CheckOrderAgrees<3>(f);
        CheckOrderAgrees<4>(f);
    }
    {
        const CubicSpline f{kNodes, y, BoundaryCondition::NotAKnot, 0.0, 0.0};
        CheckOrderAgrees<0>(f);
        CheckOrderAgrees<1>(f);
        CheckOrderAgrees<2>(f);
        CheckOrderAgrees<3>(f);
    }
    {
        const AkimaSpline f{kNodes, y};
        CheckOrderAgrees<0>(f);
        CheckOrderAgrees<1>(f);
        CheckOrderAgrees<2>(f);
        CheckOrderAgrees<3>(f);
        CheckOrderAgrees<4>(f);
    }
}

TEST(EvaluateAtNodes, ReproducesTheOrdinatesExactly) {
    const auto y = SampledOn(kNodes);
    const CubicSpline spline{kNodes, y};
    const AkimaSpline akima{kNodes, y};
    const Linear linear{kNodes, y};

    // The interpolation property, read straight off the nodal sweep.
    for (const auto &values :
         {spline.NodeValues(), akima.NodeValues(), linear.NodeValues()}) {
        ASSERT_EQ(values.size(), y.size());
        for (std::size_t k = 0; k < y.size(); ++k) {
            EXPECT_NEAR(values[k], y[k], 1.0e-13);
        }
    }
}

TEST(EvaluateAtNodes, TheTwoSidesAgreeBelowTheDiscontinuousOrder) {
    const auto y = SampledOn(kNodes);
    const CubicSpline spline{kNodes, y};

    for (std::size_t order = 0; order < 3; ++order) {
        std::vector<double> right(kNodes.size()), left(kNodes.size());
        const auto sweep = [&](Side side, std::vector<double> &out) {
            switch (order) {
            case 0:
                spline.EvaluateAtNodes<0>(std::span<double>{out}, side);
                break;
            case 1:
                spline.EvaluateAtNodes<1>(std::span<double>{out}, side);
                break;
            default:
                spline.EvaluateAtNodes<2>(std::span<double>{out}, side);
                break;
            }
        };
        sweep(Side::Right, right);
        sweep(Side::Left, left);

        for (std::size_t k = 0; k < kNodes.size(); ++k) {
            SCOPED_TRACE("order " + std::to_string(order) + ", node " +
                         std::to_string(k));
            // A cubic spline is C2, so only rounding separates the two limits
            // up to the second derivative.
            EXPECT_NEAR(left[k], right[k], 1.0e-12);
        }
    }
}

TEST(EvaluateAtNodes, TheSidesStraddleTheJumpInTheThirdDerivative) {
    const auto y = SampledOn(kNodes);
    const CubicSpline spline{kNodes, y};
    const auto n = kNodes.size();

    const auto right = spline.NodeValues<3>(Side::Right);
    const auto left = spline.NodeValues<3>(Side::Left);

    // The third derivative is constant on each segment, so its value anywhere
    // inside a segment is exactly the one-sided limit at either end. That is
    // what makes this an equality rather than a limit.
    for (std::size_t k = 0; k < n; ++k) {
        SCOPED_TRACE("node " + std::to_string(k));
        const auto rightSegment = k + 1 < n ? k : n - 2;
        const auto leftSegment = k > 0 ? k - 1 : 0;
        EXPECT_EQ(right[k], spline.Evaluate<3>(Inside(rightSegment)));
        EXPECT_EQ(left[k], spline.Evaluate<3>(Inside(leftSegment)));
    }

    // And they really do differ, so the choice of side is not cosmetic.
    bool anyJump = false;
    for (std::size_t k = 1; k + 1 < n; ++k) {
        anyJump = anyJump || left[k] != right[k];
    }
    EXPECT_TRUE(anyJump)
        << "the third derivative was continuous at every knot, so this test "
           "is not exercising the sided path";
}

TEST(EvaluateAtNodes, TheFirstDerivativeOfALinearInterpolantIsSided) {
    const auto y = SampledOn(kNodes);
    const Linear f{kNodes, y};
    const auto n = kNodes.size();

    const auto right = f.NodeValues<1>(Side::Right);
    const auto left = f.NodeValues<1>(Side::Left);

    for (std::size_t k = 0; k < n; ++k) {
        SCOPED_TRACE("node " + std::to_string(k));
        const auto rightSegment = k + 1 < n ? k : n - 2;
        const auto leftSegment = k > 0 ? k - 1 : 0;
        EXPECT_EQ(right[k], f.Evaluate<1>(Inside(rightSegment)));
        EXPECT_EQ(left[k], f.Evaluate<1>(Inside(leftSegment)));
    }
    EXPECT_NE(left[1], right[1]);
}

TEST(EvaluateAtNodes, WorksForComplexOrdinates) {
    using C = std::complex<double>;
    std::vector<C> y;
    y.reserve(kNodes.size());
    for (const auto value : kNodes) {
        y.emplace_back(std::cos(1.1 * value), std::sin(0.7 * value));
    }

    const CubicSpline spline{kNodes, y};
    const auto values = spline.NodeValues<1>();
    for (std::size_t k = 0; k < kNodes.size(); ++k) {
        EXPECT_EQ(values[k], spline.Evaluate<1>(kNodes[k]));
    }
}

TEST(EvaluateAtNodes, RejectsAnOutputOfTheWrongLength) {
    const auto y = SampledOn(kNodes);
    std::vector<double> tooShort(kNodes.size() - 1);

    const Linear linear{kNodes, y};
    const CubicSpline spline{kNodes, y};
    const AkimaSpline akima{kNodes, y};

    EXPECT_THROW(linear.EvaluateAtNodes(std::span<double>{tooShort}),
                 std::invalid_argument);
    EXPECT_THROW(spline.EvaluateAtNodes(std::span<double>{tooShort}),
                 std::invalid_argument);
    EXPECT_THROW(akima.EvaluateAtNodes(std::span<double>{tooShort}),
                 std::invalid_argument);
}

TEST(EvaluateAtNodes, TheSweepAllocatesNothing) {
    if constexpr (!InterpolationTest::AllocationCountingEnabled()) {
        GTEST_SKIP() << "allocation counting stands down under sanitizers";
    } else {
        {
            const auto before = InterpolationTest::AllocationCount();
            void *probe = ::operator new(64);
            ::operator delete(probe);
            ASSERT_GT(InterpolationTest::AllocationCount(), before)
                << "the allocation counter is not observing allocations";
        }

        const auto y = SampledOn(kNodes);
        const CubicSpline spline{kNodes, y};
        const AkimaSpline akima{kNodes, y};
        const Linear linear{kNodes, y};
        std::vector<double> out(kNodes.size());

        const auto before = InterpolationTest::AllocationCount();
        for (int i = 0; i < 100; ++i) {
            spline.EvaluateAtNodes<1>(std::span<double>{out});
            akima.EvaluateAtNodes<1>(std::span<double>{out});
            linear.EvaluateAtNodes<0>(std::span<double>{out});
        }
        const auto after = InterpolationTest::AllocationCount();

        EXPECT_EQ(after, before)
            << "the nodal sweep allocated " << (after - before) << " times";
    }
}

TEST(EvaluateAtNodes, TheOperatorShapeReproducesASplinePerLine) {
    // The whole point, end to end: one factorisation and one pair of buffers,
    // swept over many datasets on a shared grid, must give exactly what
    // building a spline per dataset gives.
    const auto seed = InterpolationTest::ReportedSeed(20260823);
    const auto n = kNodes.size();
    const CubicSplineSystem system{kNodes, BoundaryCondition::NotAKnot};

    std::vector<double> curvature(n), derivative(n);
    std::mt19937_64 engine{seed};
    std::uniform_real_distribution<double> value{-2.0, 2.0};

    for (int line = 0; line < 32; ++line) {
        std::vector<double> y(n);
        for (auto &entry : y) {
            entry = value(engine);
        }

        system.Solve(y, std::span<double>{curvature});
        system.EvaluateAtNodes<1>(y, curvature, std::span<double>{derivative});

        const CubicSpline reference{kNodes, y, BoundaryCondition::NotAKnot, 0.0,
                                    0.0};
        for (std::size_t k = 0; k < n; ++k) {
            SCOPED_TRACE("line " + std::to_string(line) + ", node " +
                         std::to_string(k));
            EXPECT_EQ(derivative[k], reference.Evaluate<1>(kNodes[k]));
        }
    }
}
