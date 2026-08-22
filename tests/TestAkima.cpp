#include <gtest/gtest.h>

#include <Interpolation/AkimaSpline.hpp>
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

template <typename x_value_t, typename y_value_t>
std::pair<y_value_t, y_value_t>
ReferenceHermite(const std::vector<x_value_t> &x,
                 const std::vector<y_value_t> &y,
                 const std::vector<y_value_t> &slopes, x_value_t query) {
    const auto upperIterator = std::ranges::upper_bound(x, query);
    auto upper =
        static_cast<std::size_t>(std::distance(x.begin(), upperIterator));
    upper = std::clamp(upper, std::size_t{1}, x.size() - 1);
    const auto lower = upper - 1;
    const auto h = x[upper] - x[lower];
    const auto secant = (y[upper] - y[lower]) / h;
    const auto c = (static_cast<y_value_t>(3) * secant -
                    static_cast<y_value_t>(2) * slopes[lower] - slopes[upper]) /
                   h;
    const auto d =
        (slopes[lower] + slopes[upper] - static_cast<y_value_t>(2) * secant) /
        (h * h);
    const auto offset = query - x[lower];
    const auto value =
        y[lower] + offset * (slopes[lower] + offset * (c + d * offset));
    const auto derivative =
        slopes[lower] + offset * (static_cast<y_value_t>(2) * c +
                                  static_cast<y_value_t>(3) * d * offset);
    return {value, derivative};
}

template <Interpolation::Real real_t>
void
CheckComplexLinearFunction() {
    using Complex = std::complex<real_t>;
    const std::vector<real_t> x{-1, 0, 2, 5, 9};
    const Complex intercept{1, 2};
    const Complex slope{2, static_cast<real_t>(-0.5)};
    std::vector<Complex> y;
    for (const auto value : x) {
        y.push_back(intercept + slope * value);
    }
    const Interpolation::AkimaSpline akima{x, y};

    for (const real_t query :
         {static_cast<real_t>(-2), static_cast<real_t>(-1),
          static_cast<real_t>(0.5), static_cast<real_t>(4),
          static_cast<real_t>(9), static_cast<real_t>(10)}) {
        InterpolationTest::ExpectScaledNear(akima(query),
                                            intercept + slope * query);
        InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(query),
                                            slope);
    }
}

} // namespace

TEST(AkimaSpline, NonuniformKnownAnswerAtBoundariesAndOutsideDomain) {
    const std::vector<double> x{0.0, 1.0, 3.0, 6.0, 10.0, 15.0};
    const std::vector<double> y{0.0, 1.0, 0.0, 2.0, 1.0, 3.0};
    const std::vector<double> expectedSlopes{
        1.0, 1.0 / 4.0, 13.0 / 58.0, 17.0 / 218.0, 3.0 / 40.0, 2.0 / 5.0};
    const Interpolation::AkimaSpline akima{x, y};

    for (const double query : {-2.0, 0.0, 0.5, 1.0, 2.0, 3.0, 4.5, 6.0, 8.0,
                               10.0, 12.5, 15.0, 17.0}) {
        const auto expected = ReferenceHermite(x, y, expectedSlopes, query);
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(akima(query), expected.first);
        InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(query),
                                            expected.second);
    }

    for (std::size_t i = 0; i < x.size(); ++i) {
        InterpolationTest::ExpectScaledNear(akima(x[i]), y[i]);
    }
    InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(x.front()),
                                        expectedSlopes.front());
    InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(x.back()),
                                        expectedSlopes.back());
}

TEST(AkimaSpline, LinearFunctionIsRecovered) {
    const std::vector<double> x{-1.0, 0.0, 2.0, 5.0, 9.0};
    const double intercept = 1.0;
    const double slope = 2.0;
    std::vector<double> y;
    for (const auto value : x) {
        y.push_back(intercept + slope * value);
    }
    const Interpolation::AkimaSpline akima{x, y};

    for (const double query :
         {-2.0, -1.0, -0.5, 0.0, 1.0, 4.0, 8.5, 9.0, 10.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(akima(query),
                                            intercept + slope * query);
        InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(query),
                                            slope);
    }
}

TEST(AkimaSpline, ComplexLinearFunctionIsRecovered) {
    CheckComplexLinearFunction<float>();
    CheckComplexLinearFunction<double>();
    CheckComplexLinearFunction<long double>();
}

TEST(AkimaSpline, ComplexNonuniformKnownAnswer) {
    using Complex = std::complex<double>;
    const std::vector<double> x{0.0, 1.0, 3.0, 6.0, 10.0, 15.0};
    const std::vector<double> realY{0.0, 1.0, 0.0, 2.0, 1.0, 3.0};
    const std::vector<double> realSlopes{1.0,          1.0 / 4.0,  13.0 / 58.0,
                                         17.0 / 218.0, 3.0 / 40.0, 2.0 / 5.0};
    std::vector<Complex> y;
    std::vector<Complex> expectedSlopes;
    for (std::size_t i = 0; i < x.size(); ++i) {
        y.emplace_back(realY[i], 2.0 + 0.5 * x[i]);
        expectedSlopes.emplace_back(realSlopes[i], 0.5);
    }
    const Interpolation::AkimaSpline akima{x, y};

    for (const double query :
         {-2.0, 0.0, 0.5, 2.0, 4.5, 8.0, 12.5, 15.0, 17.0}) {
        const auto expected = ReferenceHermite(x, y, expectedSlopes, query);
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(akima(query), expected.first);
        InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(query),
                                            expected.second);
    }
}

TEST(AkimaSpline, ThreeNodeNonlinearInputUsesEndpointSlopes) {
    const std::vector<double> x{0.0, 1.0, 3.0};
    const std::vector<double> y{0.0, 1.0, 0.0};
    const std::vector<double> expectedSlopes{1.0, 0.25, -0.5};
    const Interpolation::AkimaSpline akima{x, y};

    for (const double query : {-1.0, 0.0, 0.5, 1.0, 2.5, 3.0, 4.0}) {
        const auto expected = ReferenceHermite(x, y, expectedSlopes, query);
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(akima(query), expected.first);
        InterpolationTest::ExpectScaledNear(akima.template Evaluate<1>(query),
                                            expected.second);
    }
}

TEST(AkimaSpline, TwoPointInputIsRejected) {
    const std::vector<double> x{0.0, 1.0};
    const std::vector<double> y{0.0, 1.0};
    EXPECT_THROW((Interpolation::AkimaSpline{x, y}), std::invalid_argument);
}

TEST(AkimaSpline, MismatchedLengthsAreRejected) {
    const std::vector<double> x{0.0, 1.0, 2.0};
    const std::vector<double> y{0.0, 1.0};
    EXPECT_THROW((Interpolation::AkimaSpline{x, y}), std::invalid_argument);
}

TEST(AkimaSpline, RepeatedAbscissaIsRejected) {
    // A doubled node is not an encoding for a discontinuity here: the base
    // interpolators require strictly increasing abscissae.
    const std::vector<double> x{0.0, 1.0, 1.0, 2.0};
    const std::vector<double> y{0.0, 1.0, 2.0, 3.0};
    EXPECT_THROW((Interpolation::AkimaSpline{x, y}), std::invalid_argument);
}

TEST(AkimaSpline, UnsortedAbscissaeAreRejected) {
    const std::vector<double> x{0.0, 2.0, 1.0};
    const std::vector<double> y{0.0, 1.0, 2.0};
    EXPECT_THROW((Interpolation::AkimaSpline{x, y}), std::invalid_argument);
}
