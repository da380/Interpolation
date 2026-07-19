#include <gtest/gtest.h>

#include <Interpolation/Lagrange>
#include <complex>
#include <vector>

#include "TestUtilities.h"

TEST(LagrangePolynomial, CardinalValuesPartitionAndDerivatives) {
    const std::vector<double> nodes{-1.0, 0.5, 2.0};
    const Interpolation::LagrangePolynomial basis{nodes.begin(), nodes.end()};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            SCOPED_TRACE(i);
            SCOPED_TRACE(j);
            InterpolationTest::ExpectScaledNear(basis(i, nodes[j]),
                                                i == j ? 1.0 : 0.0);
        }
    }

    for (const double query : {-2.0, -0.25, 1.0, 3.0}) {
        const auto expected0 = (query - 0.5) * (query - 2.0) / 4.5;
        const auto expected1 = -(query + 1.0) * (query - 2.0) / 2.25;
        const auto expected2 = (query + 1.0) * (query - 0.5) / 4.5;
        const auto derivative0 = (2.0 * query - 2.5) / 4.5;
        const auto derivative1 = -(2.0 * query - 1.0) / 2.25;
        const auto derivative2 = (2.0 * query + 0.5) / 4.5;

        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(basis(0, query), expected0);
        InterpolationTest::ExpectScaledNear(basis(1, query), expected1);
        InterpolationTest::ExpectScaledNear(basis(2, query), expected2);
        InterpolationTest::ExpectScaledNear(basis.Derivative(0, query),
                                            derivative0);
        InterpolationTest::ExpectScaledNear(basis.Derivative(1, query),
                                            derivative1);
        InterpolationTest::ExpectScaledNear(basis.Derivative(2, query),
                                            derivative2);
        InterpolationTest::ExpectScaledNear(
            basis(0, query) + basis(1, query) + basis(2, query), 1.0);
        InterpolationTest::ExpectScaledNear(
            basis.Derivative(0, query) + basis.Derivative(1, query) +
                basis.Derivative(2, query),
            0.0);
    }
}

TEST(LagrangePolynomial, OneNodeBasisIsConstant) {
    const std::vector<double> nodes{2.0};
    const Interpolation::LagrangePolynomial basis{nodes.begin(), nodes.end()};

    for (const double query : {-3.0, 2.0, 8.0}) {
        InterpolationTest::ExpectScaledNear(basis(0, query), 1.0);
        InterpolationTest::ExpectScaledNear(basis.Derivative(0, query), 0.0);
    }
}

TEST(Lagrange, NonuniformRealPolynomialAndDerivativeAreRecovered) {
    const std::vector<double> x{-1.0, 0.25, 2.0, 4.0};
    const auto polynomial = [](double value) {
        return 1.0 - 2.0 * value + 0.5 * value * value +
               1.25 * value * value * value;
    };
    const auto derivative = [](double value) {
        return -2.0 + value + 3.75 * value * value;
    };

    std::vector<double> y;
    for (const auto value : x) {
        y.push_back(polynomial(value));
    }
    const Interpolation::Lagrange interpolant{x.begin(), x.end(), y.begin()};

    for (const double query : {-2.0, -1.0, 0.5, 2.0, 3.0, 4.0, 5.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(interpolant(query),
                                            polynomial(query));
        InterpolationTest::ExpectScaledNear(interpolant.Derivative(query),
                                            derivative(query));
    }
}

TEST(Lagrange, ComplexPolynomialAndDerivativeAreRecovered) {
    using Complex = std::complex<double>;
    const std::vector<double> x{-1.0, 0.25, 2.0, 4.0};
    const Complex c0{1.0, 1.0};
    const Complex c1{-2.0, 0.5};
    const Complex c2{0.25, -1.0};
    const auto polynomial = [&](double value) {
        return c0 + c1 * value + c2 * value * value;
    };
    const auto derivative = [&](double value) {
        return c1 + 2.0 * c2 * value;
    };

    std::vector<Complex> y;
    for (const auto value : x) {
        y.push_back(polynomial(value));
    }
    const Interpolation::Lagrange interpolant{x.begin(), x.end(), y.begin()};

    for (const double query : {-2.0, -1.0, 0.5, 3.0, 5.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(interpolant(query),
                                            polynomial(query));
        InterpolationTest::ExpectScaledNear(interpolant.Derivative(query),
                                            derivative(query));
    }
}
