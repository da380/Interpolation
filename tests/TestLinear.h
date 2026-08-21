#ifndef INTERPOLATION_TEST_LINEAR_GUARD_H
#define INTERPOLATION_TEST_LINEAR_GUARD_H

#include <Interpolation/Interpolation.hpp>
#include <complex>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

template <Interpolation::RealFloatingPoint x_value_t,
          Interpolation::RealOrComplexFloatingPoint y_value_t>
int LinearCheck() {
  using namespace Interpolation;

  // Make a random linear polynomial.
  auto p = Polynomial1D<y_value_t>::Random(1);

  // set the number of sampling points randomly
  std::random_device rd{};
  std::mt19937_64 gen{rd()};
  std::uniform_int_distribution dint{5, 100};
  auto nSample = dint(gen);

  // Set the x values
  const x_value_t x1 = 0;
  const x_value_t x2 = 1;
  auto h = (x2 - x1) / static_cast<x_value_t>(nSample - 1);
  std::vector<x_value_t> x;
  std::generate_n(std::back_inserter(x), nSample,
                  [x1, h, m = 0]() mutable { return x1 + h * m++; });

  // add random shifts so that points are not equally spaced
  std::uniform_real_distribution<x_value_t> hDist(-0.1 * h, 0.1 * h);
  std::transform(std::next(x.begin()), std::prev(x.end()), std::next(x.begin()),
                 [&](auto x) { return x + hDist(gen); });

  // Set the y-values
  std::vector<y_value_t> y;
  std::transform(x.begin(), x.end(), std::back_inserter(y),
                 [&](auto x) { return p(x); });

  // Form the interpolating function.
  auto f = Linear(x.begin(), x.end(), y.begin());

  // Compare exact and interpolated values at randomly sampled points
  std::uniform_real_distribution<x_value_t> xDist{x1, x2};
  constexpr auto eps = 1000 * std::numeric_limits<x_value_t>::epsilon();
  const int nRandom = 100;
  int count = 0;
  while (count++ < nRandom) {
    auto xx = xDist(gen);
    x_value_t functionError = std::abs(f(xx) - p(xx));
    x_value_t derivativeError = std::abs(f.Derivative(xx) - p.Derivative(xx));
    if (functionError > eps) return 1;
    if (derivativeError > eps) return 1;
  }

  return 0;
}

#endif  // INTERPOLATION_TEST_LINEAR_GUARD_H
