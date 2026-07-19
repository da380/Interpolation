#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>
#include <vector>

#include "TestCubicSpline.h"
#include "TestLinear.h"

namespace {

template <typename value_t>
void ExpectScaledNear(const value_t &actual, const value_t &expected) {
  using real_t = decltype(std::abs(expected));
  const auto scale = std::max<real_t>({1, std::abs(actual), std::abs(expected)});
  const auto tolerance = 1000 * std::numeric_limits<real_t>::epsilon() * scale;
  EXPECT_LE(std::abs(actual - expected), tolerance);
}

template <typename value_t>
void CheckAgainstDenseReference(const std::vector<double> &x,
                                const std::vector<value_t> &y, value_t ypl,
                                value_t ypr) {
  using Interpolation::CubicSpline;
  using Interpolation::CubicSplineBC;
  using CubicSplineTest::DenseReferenceSpline;
  using CubicSplineTest::RecoverIntervalSecondDerivatives;

  constexpr std::array boundaryConditions{
      std::pair{CubicSplineBC::Free, CubicSplineBC::Free},
      std::pair{CubicSplineBC::Free, CubicSplineBC::Clamped},
      std::pair{CubicSplineBC::Clamped, CubicSplineBC::Free},
      std::pair{CubicSplineBC::Clamped, CubicSplineBC::Clamped}};

  for (const auto &[left, right] : boundaryConditions) {
    SCOPED_TRACE(static_cast<int>(left));
    SCOPED_TRACE(static_cast<int>(right));
    auto spline = CubicSpline(x.begin(), x.end(), y.begin(), left, ypl, right, ypr);
    const DenseReferenceSpline reference{x, y, left, ypl, right, ypr};

    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
      for (double fraction : {0.0, 0.25, 0.5, 0.75}) {
        const auto value = x[i] + fraction * (x[i + 1] - x[i]);
        ExpectScaledNear(spline(value), reference(value));
        ExpectScaledNear(spline.Derivative(value), reference.Derivative(value));
      }
    }
    ExpectScaledNear(spline(x.back()), reference(x.back()));
    ExpectScaledNear(spline.Derivative(x.back()), reference.Derivative(x.back()));

    const auto first = RecoverIntervalSecondDerivatives(
        spline, x[0], x[1], y[0], y[1]);
    const auto lastIndex = x.size() - 1;
    const auto last = RecoverIntervalSecondDerivatives(
        spline, x[lastIndex - 1], x[lastIndex], y[lastIndex - 1], y[lastIndex]);
    if (left == CubicSplineBC::Free) {
      ExpectScaledNear(first.first, value_t{});
    } else {
      ExpectScaledNear(spline.Derivative(x.front()), ypl);
    }
    if (right == CubicSplineBC::Free) {
      ExpectScaledNear(last.second, value_t{});
    } else {
      ExpectScaledNear(spline.Derivative(x.back()), ypr);
    }
  }
}

}  // namespace

// Tests for linear interpolation

TEST(Linear, CheckRealSingle) {
  int i = LinearCheck<float, float>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckRealDouble) {
  int i = LinearCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckRealLongDouble) {
  int i = LinearCheck<long double, long double>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexSingle) {
  int i = LinearCheck<float, std::complex<float>>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexDouble) {
  int i = LinearCheck<double, std::complex<double>>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexLongDouble) {
  int i = LinearCheck<long double, std::complex<long double>>();
  EXPECT_EQ(0, i);
}

// Tests for cubic spline interpolation

TEST(CubicSpline, CheckRealSingle) {
  int i = CubicSplineCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckRealDouble) {
  int i = CubicSplineCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckRealLongDouble) {
  int i = CubicSplineCheck<long double, long double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexSingle) {
  int i = CubicSplineCheck<float, std::complex<float>>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexDouble) {
  int i = CubicSplineCheck<double, std::complex<double>>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexLongDouble) {
  int i = CubicSplineCheck<long double, std::complex<long double>>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, NaturalNonuniformKnownAnswer) {
  const std::vector<double> x{0.0, 1.0, 3.0, 4.0};
  const std::vector<double> y{0.0, 1.0, 0.0, 2.0};
  const Interpolation::CubicSpline spline{x.begin(), x.end(), y.begin()};

  struct Sample {
    double x;
    double value;
    double derivative;
  };
  constexpr std::array samples{
      Sample{0.5, 0.6640625, 1.109375},
      Sample{2.0, 0.3125, -1.0},
      Sample{3.5, 0.7890625, 2.140625}};

  for (const auto &sample : samples) {
    ExpectScaledNear(spline(sample.x), sample.value);
    ExpectScaledNear(spline.Derivative(sample.x), sample.derivative);
  }

  const auto first = CubicSplineTest::RecoverIntervalSecondDerivatives(
      spline, x[0], x[1], y[0], y[1]);
  const auto last = CubicSplineTest::RecoverIntervalSecondDerivatives(
      spline, x[2], x[3], y[2], y[3]);
  ExpectScaledNear(first.first, 0.0);
  ExpectScaledNear(last.second, 0.0);
}

TEST(CubicSpline, BoundaryCombinationsMatchDenseReference) {
  CheckAgainstDenseReference<double>({0.0, 1.0, 3.0, 4.0},
                                     {0.0, 1.0, 0.0, 2.0}, -0.5, 1.25);
}

TEST(CubicSpline, TwoPointBoundaryCombinationsMatchDenseReference) {
  CheckAgainstDenseReference<double>({0.0, 2.0}, {1.0, 4.0}, -0.5, 1.25);
}

TEST(CubicSpline, ComplexBoundaryCombinationsMatchDenseReference) {
  using Complex = std::complex<double>;
  CheckAgainstDenseReference<Complex>(
      {0.0, 1.0, 3.0, 4.0},
      {Complex{0.0, 1.0}, Complex{1.0, -0.5}, Complex{0.0, 0.25},
       Complex{2.0, -1.0}},
      Complex{-0.5, 0.75}, Complex{1.25, -0.25});
}
