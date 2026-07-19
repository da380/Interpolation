#include <gtest/gtest.h>

#include <Interpolation/Polynomial>
#include <complex>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

template <typename value_t>
void ExpectCoefficients(
    const Interpolation::Polynomial1D<value_t> &polynomial,
    const std::vector<value_t> &expected) {
    ASSERT_EQ(polynomial.polycoeff().size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        SCOPED_TRACE(i);
        InterpolationTest::ExpectScaledNear(polynomial.polycoeff(i),
                                            expected[i]);
    }
}

}   // namespace

TEST(Polynomial1D, DefaultConstructionIsZeroPolynomial) {
    Interpolation::Polynomial1D<double> zero;

    EXPECT_EQ(zero.Degree(), 0);
    ExpectCoefficients(zero, std::vector<double>{0.0});
    InterpolationTest::ExpectScaledNear(zero(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.Derivative(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.Primitive(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.Integrate(-1.0, 2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero[0], 0.0);
    InterpolationTest::ExpectScaledNear(zero.polycoeff(-1), 0.0);
    InterpolationTest::ExpectScaledNear(zero.polycoeff(1), 0.0);

    const Interpolation::Polynomial1D<double> polynomial{1.0, -2.0, 3.0};
    ExpectCoefficients(zero + polynomial,
                       std::vector<double>{1.0, -2.0, 3.0});
    ExpectCoefficients(zero * polynomial,
                       std::vector<double>{0.0, 0.0, 0.0});
    ExpectCoefficients(zero + 2.0, std::vector<double>{2.0});

    const Interpolation::Polynomial1D<std::complex<double>> complexZero;
    ExpectCoefficients(complexZero,
                       std::vector<std::complex<double>>{{0.0, 0.0}});
}

TEST(Polynomial1D, ConstructsFromConstCoefficientVector) {
    const std::vector<double> coefficients{1.0, -2.0, 3.0};
    const Interpolation::Polynomial1D<double> polynomial{coefficients};

    ExpectCoefficients(polynomial, coefficients);
}

TEST(Polynomial1D, EvaluationCalculusAndCoefficientAccess) {
    Interpolation::Polynomial1D<double> polynomial{1.0, -2.0, 3.0};

    EXPECT_EQ(polynomial.Degree(), 2);
    InterpolationTest::ExpectScaledNear(polynomial(2.0), 9.0);
    InterpolationTest::ExpectScaledNear(polynomial.Derivative(2.0), 10.0);
    InterpolationTest::ExpectScaledNear(polynomial.Primitive(2.0), 6.0);
    InterpolationTest::ExpectScaledNear(polynomial.Integrate(-1.0, 2.0), 9.0);

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    InterpolationTest::ExpectScaledNear(polynomial.Primative(2.0), 6.0);
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

    InterpolationTest::ExpectScaledNear(polynomial[1], -2.0);
    InterpolationTest::ExpectScaledNear(polynomial.polycoeff(-1), 0.0);
    InterpolationTest::ExpectScaledNear(polynomial.polycoeff(3), 0.0);
    ExpectCoefficients(polynomial, std::vector<double>{1.0, -2.0, 3.0});

    *polynomial.begin() = 4.0;
    InterpolationTest::ExpectScaledNear(polynomial(0.0), 4.0);
    EXPECT_EQ(std::distance(polynomial.cbegin(), polynomial.cend()), 3);
}

TEST(Polynomial1D, ScalarAndPolynomialArithmetic) {
    const Interpolation::Polynomial1D<double> polynomial{1.0, -2.0, 3.0};
    const Interpolation::Polynomial1D<double> other{-1.0, 4.0};

    ExpectCoefficients(polynomial + 2.0,
                       std::vector<double>{3.0, -2.0, 3.0});
    ExpectCoefficients(2.0 + polynomial,
                       std::vector<double>{3.0, -2.0, 3.0});
    ExpectCoefficients(polynomial - 2.0,
                       std::vector<double>{-1.0, -2.0, 3.0});
    ExpectCoefficients(2.0 - polynomial,
                       std::vector<double>{1.0, 2.0, -3.0});
    ExpectCoefficients(polynomial * 2.0,
                       std::vector<double>{2.0, -4.0, 6.0});
    ExpectCoefficients(2.0 * polynomial,
                       std::vector<double>{2.0, -4.0, 6.0});
    ExpectCoefficients(polynomial / 2.0,
                       std::vector<double>{0.5, -1.0, 1.5});

    ExpectCoefficients(polynomial + other,
                       std::vector<double>{0.0, 2.0, 3.0});
    ExpectCoefficients(polynomial - other,
                       std::vector<double>{2.0, -6.0, 3.0});
    ExpectCoefficients(polynomial * other,
                       std::vector<double>{-1.0, 6.0, -11.0, 12.0});
    ExpectCoefficients(-polynomial,
                       std::vector<double>{-1.0, 2.0, -3.0});

    auto compoundScalar = polynomial;
    compoundScalar += 2.0;
    compoundScalar -= 1.0;
    compoundScalar *= 4.0;
    compoundScalar /= 2.0;
    ExpectCoefficients(compoundScalar,
                       std::vector<double>{4.0, -4.0, 6.0});

    auto compoundPolynomial = polynomial;
    compoundPolynomial += other;
    compoundPolynomial -= other;
    compoundPolynomial *= other;
    ExpectCoefficients(compoundPolynomial,
                       std::vector<double>{-1.0, 6.0, -11.0, 12.0});
}

TEST(Polynomial1D, ConversionComplexValuesAndStreamOutput) {
    Interpolation::Polynomial1D<float> source{1.0F, 2.0F, -0.5F};
    const Interpolation::Polynomial1D<double> converted{source};
    ExpectCoefficients(converted, std::vector<double>{1.0, 2.0, -0.5});

    using Complex = std::complex<double>;
    const Interpolation::Polynomial1D<Complex> complexPolynomial{
        Complex{1.0, 2.0}, Complex{2.0, -1.0}};
    const Complex query{2.0, 0.0};
    InterpolationTest::ExpectScaledNear(complexPolynomial(query),
                                        Complex{5.0, 0.0});
    InterpolationTest::ExpectScaledNear(complexPolynomial.Derivative(query),
                                        Complex{2.0, -1.0});
    InterpolationTest::ExpectScaledNear(complexPolynomial.Primitive(query),
                                        Complex{6.0, 2.0});
    InterpolationTest::ExpectScaledNear(
        complexPolynomial.Integrate(Complex{0.0, 0.0}, query),
        Complex{6.0, 2.0});

    std::ostringstream output;
    output << converted;
    EXPECT_EQ(output.str(), "1 2 -0.5 ");
}

TEST(Polynomial1D, CopyMoveAndCrossTypeAssignment) {
    using FloatPolynomial = Interpolation::Polynomial1D<float>;
    using DoublePolynomial = Interpolation::Polynomial1D<double>;
    using ComplexPolynomial =
        Interpolation::Polynomial1D<std::complex<double>>;
    static_assert(std::is_copy_assignable_v<DoublePolynomial>);
    static_assert(std::is_move_assignable_v<DoublePolynomial>);

    const FloatPolynomial source{1.0F, 2.0F, -0.5F};
    DoublePolynomial converted{9.0, 8.0};
    converted = source;
    ExpectCoefficients(converted, std::vector<double>{1.0, 2.0, -0.5});
    ExpectCoefficients(source, std::vector<float>{1.0F, 2.0F, -0.5F});

    DoublePolynomial copied{7.0};
    copied = converted;
    ExpectCoefficients(copied, std::vector<double>{1.0, 2.0, -0.5});

    DoublePolynomial moveSource{4.0, -3.0, 2.0, -1.0};
    DoublePolynomial moved{0.0};
    moved = std::move(moveSource);
    ExpectCoefficients(moved,
                       std::vector<double>{4.0, -3.0, 2.0, -1.0});

    ComplexPolynomial complexAssigned{{0.0, 1.0}};
    complexAssigned = converted;
    ExpectCoefficients(complexAssigned,
                       std::vector<std::complex<double>>{{1.0, 0.0},
                                                         {2.0, 0.0},
                                                         {-0.5, 0.0}});
}
