#include <gtest/gtest.h>

#include <Interpolation/Polynomial.hpp>
#include <complex>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

#include "TestUtilities.h"

namespace {

template <typename value_t>
void
ExpectCoefficients(const Interpolation::Polynomial<value_t> &polynomial,
                   const std::vector<value_t> &expected) {
    ASSERT_EQ(polynomial.polycoeff().size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        SCOPED_TRACE(i);
        InterpolationTest::ExpectScaledNear(polynomial.polycoeff(i),
                                            expected[i]);
    }
}

} // namespace

TEST(Polynomial, DefaultConstructionIsZeroPolynomial) {
    Interpolation::Polynomial<double> zero;

    EXPECT_EQ(zero.Degree(), 0);
    ExpectCoefficients(zero, std::vector<double>{0.0});
    InterpolationTest::ExpectScaledNear(zero(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.template Evaluate<1>(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.Primitive(2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero.Integrate(-1.0, 2.0), 0.0);
    InterpolationTest::ExpectScaledNear(zero[0], 0.0);
    InterpolationTest::ExpectScaledNear(zero.polycoeff(-1), 0.0);
    InterpolationTest::ExpectScaledNear(zero.polycoeff(1), 0.0);

    const Interpolation::Polynomial<double> polynomial{1.0, -2.0, 3.0};
    ExpectCoefficients(zero + polynomial, std::vector<double>{1.0, -2.0, 3.0});
    ExpectCoefficients(zero * polynomial, std::vector<double>{0.0, 0.0, 0.0});
    ExpectCoefficients(zero + 2.0, std::vector<double>{2.0});

    const Interpolation::Polynomial<std::complex<double>> complexZero;
    ExpectCoefficients(complexZero,
                       std::vector<std::complex<double>>{{0.0, 0.0}});
}

TEST(Polynomial, ConstructsFromConstCoefficientVector) {
    const std::vector<double> coefficients{1.0, -2.0, 3.0};
    const Interpolation::Polynomial<double> polynomial{coefficients};

    ExpectCoefficients(polynomial, coefficients);
}

TEST(Polynomial, EvaluationCalculusAndCoefficientAccess) {
    Interpolation::Polynomial<double> polynomial{1.0, -2.0, 3.0};

    EXPECT_EQ(polynomial.Degree(), 2);
    InterpolationTest::ExpectScaledNear(polynomial(2.0), 9.0);
    InterpolationTest::ExpectScaledNear(polynomial.template Evaluate<1>(2.0),
                                        10.0);
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

TEST(Polynomial, ScalarAndPolynomialArithmetic) {
    const Interpolation::Polynomial<double> polynomial{1.0, -2.0, 3.0};
    const Interpolation::Polynomial<double> other{-1.0, 4.0};

    ExpectCoefficients(polynomial + 2.0, std::vector<double>{3.0, -2.0, 3.0});
    ExpectCoefficients(2.0 + polynomial, std::vector<double>{3.0, -2.0, 3.0});
    ExpectCoefficients(polynomial - 2.0, std::vector<double>{-1.0, -2.0, 3.0});
    ExpectCoefficients(2.0 - polynomial, std::vector<double>{1.0, 2.0, -3.0});
    ExpectCoefficients(polynomial * 2.0, std::vector<double>{2.0, -4.0, 6.0});
    ExpectCoefficients(2.0 * polynomial, std::vector<double>{2.0, -4.0, 6.0});
    ExpectCoefficients(polynomial / 2.0, std::vector<double>{0.5, -1.0, 1.5});

    ExpectCoefficients(polynomial + other, std::vector<double>{0.0, 2.0, 3.0});
    ExpectCoefficients(polynomial - other, std::vector<double>{2.0, -6.0, 3.0});
    ExpectCoefficients(polynomial * other,
                       std::vector<double>{-1.0, 6.0, -11.0, 12.0});
    ExpectCoefficients(-polynomial, std::vector<double>{-1.0, 2.0, -3.0});

    auto compoundScalar = polynomial;
    compoundScalar += 2.0;
    compoundScalar -= 1.0;
    compoundScalar *= 4.0;
    compoundScalar /= 2.0;
    ExpectCoefficients(compoundScalar, std::vector<double>{4.0, -4.0, 6.0});

    auto compoundPolynomial = polynomial;
    compoundPolynomial += other;
    compoundPolynomial -= other;
    compoundPolynomial *= other;
    ExpectCoefficients(compoundPolynomial,
                       std::vector<double>{-1.0, 6.0, -11.0, 12.0});
}

TEST(Polynomial, ConversionComplexValuesAndStreamOutput) {
    Interpolation::Polynomial<float> source{1.0F, 2.0F, -0.5F};
    const Interpolation::Polynomial<double> converted{source};
    ExpectCoefficients(converted, std::vector<double>{1.0, 2.0, -0.5});

    using Complex = std::complex<double>;
    const Interpolation::Polynomial<Complex> complexPolynomial{
        Complex{1.0, 2.0}, Complex{2.0, -1.0}};
    const Complex query{2.0, 0.0};
    InterpolationTest::ExpectScaledNear(complexPolynomial(query),
                                        Complex{5.0, 0.0});
    InterpolationTest::ExpectScaledNear(
        complexPolynomial.template Evaluate<1>(query), Complex{2.0, -1.0});
    InterpolationTest::ExpectScaledNear(complexPolynomial.Primitive(query),
                                        Complex{6.0, 2.0});
    InterpolationTest::ExpectScaledNear(
        complexPolynomial.Integrate(Complex{0.0, 0.0}, query),
        Complex{6.0, 2.0});

    std::ostringstream output;
    output << converted;
    EXPECT_EQ(output.str(), "1 2 -0.5 ");
}

TEST(Polynomial, CopyMoveAndCrossTypeAssignment) {
    using FloatPolynomial = Interpolation::Polynomial<float>;
    using DoublePolynomial = Interpolation::Polynomial<double>;
    using ComplexPolynomial = Interpolation::Polynomial<std::complex<double>>;
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
    ExpectCoefficients(moved, std::vector<double>{4.0, -3.0, 2.0, -1.0});

    ComplexPolynomial complexAssigned{{0.0, 1.0}};
    complexAssigned = converted;
    ExpectCoefficients(
        complexAssigned,
        std::vector<std::complex<double>>{{1.0, 0.0}, {2.0, 0.0}, {-0.5, 0.0}});
}

TEST(Polynomial, DegreeIsHonestAfterCancellation) {
    // Before trimming, a subtraction that cancels the leading terms left the
    // polynomial claiming a degree it no longer had.
    const Interpolation::Polynomial<double> a{1.0, 2.0, 3.0};
    const Interpolation::Polynomial<double> b{0.0, 0.0, 3.0};

    const auto difference = a - b;
    EXPECT_EQ(difference.Degree(), 1);
    InterpolationTest::ExpectScaledNear(difference(2.0), 5.0);

    // Cancelling everything leaves the zero polynomial, degree 0, not an
    // empty coefficient vector.
    const auto zero = a - a;
    EXPECT_EQ(zero.Degree(), 0);
    InterpolationTest::ExpectScaledNear(zero(3.0), 0.0);
}

TEST(Polynomial, MixedArithmeticPromotesRatherThanTruncates) {
    const Interpolation::Polynomial<double> real{1.0, 2.0};
    const std::complex<double> scalar{0.0, 1.0};

    // Taking the left operand's type would have silently dropped the
    // imaginary part; common_type_t promotes instead.
    const auto sum = real + scalar;
    static_assert(
        std::is_same_v<decltype(sum),
                       const Interpolation::Polynomial<std::complex<double>>>);
    EXPECT_DOUBLE_EQ(sum(0.0).real(), 1.0);
    EXPECT_DOUBLE_EQ(sum(0.0).imag(), 1.0);

    const auto scaled = scalar * real;
    static_assert(
        std::is_same_v<decltype(scaled),
                       const Interpolation::Polynomial<std::complex<double>>>);
    EXPECT_DOUBLE_EQ(scaled(1.0).imag(), 3.0);

    // Integral scalars still work and stay real.
    const auto shifted = real + 2;
    static_assert(std::is_same_v<decltype(shifted),
                                 const Interpolation::Polynomial<double>>);
    InterpolationTest::ExpectScaledNear(shifted(0.0), 3.0);
}

TEST(Polynomial, ArithmeticIsFoundByArgumentDependentLookupOnly) {
    // The operators are hidden friends, so they are reachable through ADL on
    // a Polynomial operand but are not templates sitting in the global
    // namespace taking part in unrelated overload resolution.
    const Interpolation::Polynomial<double> p{1.0, 1.0};
    const auto q = p * p;
    InterpolationTest::ExpectScaledNear(q(2.0), 9.0);
    EXPECT_EQ(q.Degree(), 2);
}

TEST(Polynomial, RandomIsReproducibleFromItsSeed) {
    const auto a = Interpolation::Polynomial<double>::Random(4, 20260821ull);
    const auto b = Interpolation::Polynomial<double>::Random(4, 20260821ull);
    const auto c = Interpolation::Polynomial<double>::Random(4, 20260822ull);

    EXPECT_EQ(a.Degree(), 4);
    for (int i = 0; i <= a.Degree(); ++i) {
        EXPECT_DOUBLE_EQ(a.polycoeff(i), b.polycoeff(i));
    }
    // A different seed should give different coefficients.
    bool differs = false;
    for (int i = 0; i <= a.Degree(); ++i) {
        differs = differs || a.polycoeff(i) != c.polycoeff(i);
    }
    EXPECT_TRUE(differs);
}

TEST(Polynomial, EvaluateMatchesRepeatedDifferentiation) {
    // p(x) = 1 + 2x + 3x^2 + 4x^3
    const Interpolation::Polynomial<double> p{1.0, 2.0, 3.0, 4.0};
    InterpolationTest::ExpectScaledNear(p.Evaluate<0>(2.0), p(2.0));
    InterpolationTest::ExpectScaledNear(p.Evaluate<1>(2.0),
                                        2.0 + 6.0 * 2.0 + 12.0 * 4.0);
    InterpolationTest::ExpectScaledNear(p.Evaluate<2>(2.0), 6.0 + 24.0 * 2.0);
    InterpolationTest::ExpectScaledNear(p.Evaluate<3>(2.0), 24.0);
    InterpolationTest::ExpectScaledNear(p.Evaluate<4>(2.0), 0.0);
}
