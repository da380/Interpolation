#include <gtest/gtest.h>

#include <Interpolation/Lagrange.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <vector>

#include "TestUtilities.h"

namespace {

// Two unequal node sets away from zero: the five Gauss-Lobatto-Legendre
// nodes of [0.35, 0.5], the spectral element in which the derivative beside
// a node was found wrong, and seven Chebyshev nodes of [1, 3].
std::vector<std::vector<double>>
UnequalNodeSets() {
    const auto a = 0.35;
    const auto b = 0.5;
    const auto s = std::sqrt(3.0 / 7.0);
    std::vector<std::vector<double>> sets{
        {a, a + (b - a) * (1.0 - s) / 2.0, (a + b) / 2.0,
         a + (b - a) * (1.0 + s) / 2.0, b},
        {}};
    const auto n = 7;
    for (int j = 0; j < n; ++j) {
        sets[1].push_back(2.0 - std::cos((2 * j + 1) * std::numbers::pi /
                                         (2.0 * n)));
    }
    return sets;
}

// Queries 1, 2, 4, 64 and 1e4 ulps either side of every node.
std::vector<double>
QueriesBesideNodes(const std::vector<double> &nodes) {
    std::vector<double> queries;
    for (const auto node : nodes) {
        for (const auto toward : {-std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity()}) {
            const auto ulp = std::nextafter(node, toward) - node;
            for (const double ulps : {1.0, 2.0, 4.0, 64.0, 1e4}) {
                queries.push_back(node + ulps * ulp);
            }
        }
    }
    return queries;
}

struct ProductForm {
    long double slope;
    long double scale;
};

// L_i'(t) from the product form, which never divides by t - x_j, in long
// double. The scale is the sum of the sizes of its terms, which is |L_i'(t)|
// unless they cancel, as at the turning point of a symmetric basis function.
ProductForm
CardinalSlope(const std::vector<double> &nodes, std::size_t i, long double t) {
    auto denominator = 1.0L;
    for (std::size_t m = 0; m < nodes.size(); ++m) {
        if (m != i) {
            denominator *= static_cast<long double>(nodes[i]) - nodes[m];
        }
    }
    auto slope = 0.0L;
    auto scale = 0.0L;
    for (std::size_t j = 0; j < nodes.size(); ++j) {
        if (j == i) {
            continue;
        }
        auto product = 1.0L;
        for (std::size_t m = 0; m < nodes.size(); ++m) {
            if (m != i && m != j) {
                product *= t - nodes[m];
            }
        }
        slope += product;
        scale += std::abs(product);
    }
    return {slope / denominator, scale / std::abs(denominator)};
}

template <typename Scalar>
void
ExpectSlopeAccurateBesideNodes(const std::vector<double> &x,
                               const std::vector<Scalar> &y) {
    const Interpolation::Lagrange interpolant{x, y};
    for (const auto query : QueriesBesideNodes(x)) {
        SCOPED_TRACE(query);
        // The sum of y_i L_i' is unchanged by a shift of the ordinates, and
        // shifting by y_0 makes the reference sum cancel less.
        std::complex<long double> expected{};
        for (std::size_t i = 0; i < x.size(); ++i) {
            const auto difference =
                std::complex<double>(y[i]) - std::complex<double>(y[0]);
            expected += std::complex<long double>(difference) *
                        CardinalSlope(x, i, query).slope;
        }
        const auto actual = std::complex<long double>(std::complex<double>(
            interpolant.template Evaluate<1>(query)));
        EXPECT_LE(std::abs(actual - expected), 1e-13L * std::abs(expected));
    }
}

} // namespace

TEST(LagrangeBasis, CardinalValuesPartitionAndDerivatives) {
    const std::vector<double> nodes{-1.0, 0.5, 2.0};
    const Interpolation::LagrangeBasis basis{nodes};

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
        InterpolationTest::ExpectScaledNear(
            basis.template Evaluate<1>(0, query), derivative0);
        InterpolationTest::ExpectScaledNear(
            basis.template Evaluate<1>(1, query), derivative1);
        InterpolationTest::ExpectScaledNear(
            basis.template Evaluate<1>(2, query), derivative2);
        InterpolationTest::ExpectScaledNear(
            basis(0, query) + basis(1, query) + basis(2, query), 1.0);
        InterpolationTest::ExpectScaledNear(
            basis.template Evaluate<1>(0, query) +
                basis.template Evaluate<1>(1, query) +
                basis.template Evaluate<1>(2, query),
            0.0);
    }
}

TEST(LagrangeBasis, OneNodeBasisIsConstant) {
    const std::vector<double> nodes{2.0};
    const Interpolation::LagrangeBasis basis{nodes};

    for (const double query : {-3.0, 2.0, 8.0}) {
        InterpolationTest::ExpectScaledNear(basis(0, query), 1.0);
        InterpolationTest::ExpectScaledNear(
            basis.template Evaluate<1>(0, query), 0.0);
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
    const Interpolation::Lagrange interpolant{x, y};

    for (const double query : {-2.0, -1.0, 0.5, 2.0, 3.0, 4.0, 5.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(interpolant(query),
                                            polynomial(query));
        InterpolationTest::ExpectScaledNear(
            interpolant.template Evaluate<1>(query), derivative(query));
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
    const auto derivative = [&](double value) { return c1 + 2.0 * c2 * value; };

    std::vector<Complex> y;
    for (const auto value : x) {
        y.push_back(polynomial(value));
    }
    const Interpolation::Lagrange interpolant{x, y};

    for (const double query : {-2.0, -1.0, 0.5, 3.0, 5.0}) {
        SCOPED_TRACE(query);
        InterpolationTest::ExpectScaledNear(interpolant(query),
                                            polynomial(query));
        InterpolationTest::ExpectScaledNear(
            interpolant.template Evaluate<1>(query), derivative(query));
    }
}

// Beside a node, the barycentric derivative divides L_i(x) - 1 by the small
// x - x_i, so it once kept only the rounding of L_i(x): 1e-2 relative error one
// ulp from a node.
TEST(LagrangeBasis, DerivativeIsAccurateBesideANode) {
    for (const auto &nodes : UnequalNodeSets()) {
        const Interpolation::LagrangeBasis basis{nodes};
        for (const auto query : QueriesBesideNodes(nodes)) {
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                SCOPED_TRACE(query);
                SCOPED_TRACE(i);
                const auto expected = CardinalSlope(nodes, i, query);
                const long double actual =
                    basis.template Evaluate<1>(i, query);
                EXPECT_LE(std::abs(actual - expected.slope),
                          1e-13L * expected.scale);
            }
        }
    }
}

// The same difference quotient, p(x) - y_k over x - x_k, goes wrong for every
// set of ordinates: 4e-1 relative error one ulp from a node.
TEST(Lagrange, DerivativeIsAccurateBesideANode) {
    for (const auto &x : UnequalNodeSets()) {
        std::vector<double> real;
        std::vector<std::complex<double>> complex;
        for (const auto value : x) {
            real.push_back(std::exp(value));
            complex.push_back(std::exp(std::complex<double>{1.0, 2.0} * value));
        }
        ExpectSlopeAccurateBesideNodes(x, real);
        ExpectSlopeAccurateBesideNodes(x, complex);
    }
}
