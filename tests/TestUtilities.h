#ifndef INTERPOLATION_TEST_UTILITIES_GUARD_H
#define INTERPOLATION_TEST_UTILITIES_GUARD_H

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace InterpolationTest {

/**
 * @brief Report the seed a randomised check will use, and return it.
 *
 * Seeding from std::random_device makes a CI failure unreproducible: the run
 * that failed used a value nobody recorded. These checks take a fixed seed
 * instead, and print it into the GoogleTest trace so a failure names the
 * input that produced it.
 *
 * Override at runtime with INTERPOLATION_TEST_SEED to re-run the same checks
 * on different data.
 */
inline std::uint64_t
ReportedSeed(std::uint64_t fixedSeed) {
    auto seed = fixedSeed;
    if (const char *override = std::getenv("INTERPOLATION_TEST_SEED")) {
        seed = std::strtoull(override, nullptr, 0);
    }
    std::cout << "[    SEED  ] " << seed << std::endl;
    return seed;
}

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
