#include <gtest/gtest.h>

#include <Interpolation/Interpolation.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

// Builds an interpolator from containers that go out of scope on return.
// Under the previous iterator-based API this returned three dangling
// iterators; taking rvalue ranges by value makes it well defined.
template <typename Interpolator>
Interpolator
MakeOwning() {
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{0.0, 1.0, 4.0, 9.0};
    return Interpolator{std::move(x), std::move(y)};
}

} // namespace

TEST(RangeSemantics, LvalueRangesAreBorrowedAndTrackMutation) {
    std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> y{0.0, 1.0, 2.0};
    const Interpolation::Linear linear{x, y};

    InterpolationTest::ExpectScaledNear(linear(0.5), 0.5);

    // A borrowing interpolator observes the container it was given. Only the
    // ordinates are changed here, so the abscissa invariant still holds.
    y[1] = 10.0;
    InterpolationTest::ExpectScaledNear(linear(0.5), 5.0);
}

TEST(RangeSemantics, RvalueRangesAreOwnedAndOutliveTheirSource) {
    const auto linear = MakeOwning<
        Interpolation::Linear<std::ranges::owning_view<std::vector<double>>,
                              std::ranges::owning_view<std::vector<double>>>>();
    InterpolationTest::ExpectScaledNear(linear(0.5), 0.5);
    EXPECT_EQ(linear.Size(), 4u);

    const auto spline = MakeOwning<Interpolation::CubicSpline<
        std::ranges::owning_view<std::vector<double>>,
        std::ranges::owning_view<std::vector<double>>>>();
    InterpolationTest::ExpectScaledNear(spline(0.0), 0.0);
    InterpolationTest::ExpectScaledNear(spline(3.0), 9.0);
}

TEST(RangeSemantics, DeductionGuidesSelectBorrowOrOwn) {
    std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> y{0.0, 1.0, 4.0};

    const Interpolation::Linear borrowing{x, y};
    const Interpolation::Linear owning{std::move(x), std::move(y)};

    static_assert(
        std::is_same_v<decltype(borrowing),
                       const Interpolation::Linear<
                           std::ranges::ref_view<std::vector<double>>,
                           std::ranges::ref_view<std::vector<double>>>>);
    static_assert(
        std::is_same_v<decltype(owning),
                       const Interpolation::Linear<
                           std::ranges::owning_view<std::vector<double>>,
                           std::ranges::owning_view<std::vector<double>>>>);
    SUCCEED();
}

TEST(RangeSemantics, HigherDerivativesAreConsistent) {
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> y{0.0, 1.0, 8.0, 27.0, 64.0};

    const Interpolation::Linear linear{x, y};
    // Piecewise linear: the second derivative and beyond vanish identically.
    InterpolationTest::ExpectScaledNear(linear.Evaluate<2>(1.5), 0.0);
    InterpolationTest::ExpectScaledNear(linear.Evaluate<7>(1.5), 0.0);

    const Interpolation::CubicSpline spline{x, y};
    // Cubic pieces: the third derivative is piecewise constant, so it must
    // agree with a difference quotient of the second derivative.
    const double h = 1.0e-6;
    const double thirdByDifference =
        (spline.Evaluate<2>(1.5 + h) - spline.Evaluate<2>(1.5 - h)) / (2 * h);
    EXPECT_NEAR(spline.Evaluate<3>(1.5), thirdByDifference, 1.0e-5);
    InterpolationTest::ExpectScaledNear(spline.Evaluate<4>(1.5), 0.0);

    const Interpolation::AkimaSpline akima{x, y};
    const double akimaThird =
        (akima.Evaluate<2>(1.5 + h) - akima.Evaluate<2>(1.5 - h)) / (2 * h);
    EXPECT_NEAR(akima.Evaluate<3>(1.5), akimaThird, 1.0e-5);
    InterpolationTest::ExpectScaledNear(akima.Evaluate<4>(1.5), 0.0);
}

TEST(RangeSemantics, EveryInterpolatorRejectsBadSamples) {
    const std::vector<double> tooShort{0.0};
    const std::vector<double> tooShortY{1.0};
    EXPECT_THROW((Interpolation::Linear{tooShort, tooShortY}),
                 std::invalid_argument);

    const std::vector<double> repeated{0.0, 1.0, 1.0, 2.0};
    const std::vector<double> repeatedY{0.0, 1.0, 2.0, 3.0};
    EXPECT_THROW((Interpolation::Linear{repeated, repeatedY}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::CubicSpline{repeated, repeatedY}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::AkimaSpline{repeated, repeatedY}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::Lagrange{repeated, repeatedY}),
                 std::invalid_argument);

    const std::vector<double> x{0.0, 1.0, 2.0};
    const std::vector<double> shortY{0.0, 1.0};
    EXPECT_THROW((Interpolation::Linear{x, shortY}), std::invalid_argument);
    EXPECT_THROW((Interpolation::CubicSpline{x, shortY}),
                 std::invalid_argument);
    EXPECT_THROW((Interpolation::Lagrange{x, shortY}), std::invalid_argument);
}
