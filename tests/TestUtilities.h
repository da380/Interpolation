#ifndef INTERPOLATION_TEST_UTILITIES_GUARD_H
#define INTERPOLATION_TEST_UTILITIES_GUARD_H

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace InterpolationTest {

template <typename value_t>
void
ExpectScaledNear(const value_t &actual, const value_t &expected) {
    using real_t = decltype(std::abs(expected));
    const auto scale =
        std::max<real_t>({real_t{1}, std::abs(actual), std::abs(expected)});
    const auto tolerance =
        real_t{1000} * std::numeric_limits<real_t>::epsilon() * scale;
    EXPECT_LE(std::abs(actual - expected), tolerance);
}

} // namespace InterpolationTest

#endif // INTERPOLATION_TEST_UTILITIES_GUARD_H
