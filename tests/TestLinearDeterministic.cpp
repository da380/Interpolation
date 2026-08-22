#include <gtest/gtest.h>

#include <Interpolation/Linear.hpp>
#include <complex>
#include <vector>

#include "TestUtilities.h"

TEST(LinearDeterministic, NonuniformValuesDerivativesAndExtrapolation) {
    const std::vector<double> x{0.0, 1.0, 3.0, 6.0};
    const std::vector<double> y{1.0, 3.0, -1.0, 5.0};
    const Interpolation::Linear linear{x, y};

    struct Sample {
        double x;
        double value;
        double derivative;
    };
    const std::vector<Sample> samples{
        {-0.5, 0.0, 2.0}, {0.0, 1.0, 2.0},  {0.5, 2.0, 2.0},
        {1.0, 3.0, -2.0}, {2.0, 1.0, -2.0}, {3.0, -1.0, 2.0},
        {4.5, 2.0, 2.0},  {6.0, 5.0, 2.0},  {7.0, 7.0, 2.0}};

    for (const auto &sample : samples) {
        SCOPED_TRACE(sample.x);
        InterpolationTest::ExpectScaledNear(linear(sample.x), sample.value);
        InterpolationTest::ExpectScaledNear(
            linear.template Evaluate<1>(sample.x), sample.derivative);
    }
}

TEST(LinearDeterministic, ComplexLinearFunctionIsRecovered) {
    using Complex = std::complex<double>;
    const std::vector<double> x{-1.0, 0.5, 2.0, 5.0};
    const Complex intercept{1.0, 2.0};
    const Complex slope{2.0, -0.5};
    std::vector<Complex> y;
    for (const auto value : x) {
        y.push_back(intercept + slope * value);
    }

    const Interpolation::Linear linear{x, y};
    for (const double query : {-2.0, -1.0, 0.0, 1.25, 5.0, 6.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(linear(query),
                                            intercept + slope * query);
        InterpolationTest::ExpectScaledNear(linear.template Evaluate<1>(query),
                                            slope);
    }
}
