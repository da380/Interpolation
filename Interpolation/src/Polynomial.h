#ifndef INTERPOLATION_POLYNOMIAL_GUARD_H
#define INTERPOLATION_POLYNOMIAL_GUARD_H

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <random>
#include <vector>

#include "Concepts.h"

namespace Interpolation {

template <typename T>
requires RealOrComplexFloatingPoint<T>
class Polynomial1D {
 public:
  // Construct from std::vector.
  Polynomial1D(std::vector<T> coef) : coef{coef} {}

  // Construct from std::initializer list.
  Polynomial1D(std::initializer_list<T> list) : coef{std::vector<T>{list}} {}

  // Return a real random polynomial of given degree.
  static Polynomial1D Random(int n) requires RealFloatingPoint<T> {
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<T> d{};
    std::vector<T> a;
    std::generate_n(std::back_inserter(a), n + 1, [&]() {
      if constexpr (RealFloatingPoint<T>) {
        return d(gen);
      }
    });
    return Polynomial1D(a);
  }

  // Return a complex random polynomial of given degree.
  static Polynomial1D Random(int n) requires ComplexFloatingPoint<T> {
    using S = typename T::value_type;
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<S> d{};
    std::vector<T> a;
    std::generate_n(std::back_inserter(a), n + 1, [&]() {
      return T{d(gen), d(gen)};
    });
    return Polynomial1D(a);
  }

  // Returns the degree.
  auto Degree() const { return coef.size() - 1; }

  // Evaluates the function.
  T operator()(T x) const {
    return std::accumulate(coef.rbegin(), coef.rend(), static_cast<T>(0),
                           [x](auto p, auto a) { return p * x + a; });
  }

  // Evaluates the derivative.
  T Derivative(T x) const {
    return std::accumulate(coef.rbegin(), std::prev(coef.rend()),
                           static_cast<T>(0),
                           [x, m = Degree()](auto p, auto a) mutable {
                             return p * x + static_cast<T>(m--) * a;
                           });
  }

  // Evaluates the primative.
  T Primative(T x) const {
    return std::accumulate(coef.rbegin(), coef.rend(), static_cast<T>(0),
                           [x, m = Degree() + 1](auto p, auto a) mutable {
                             return p * x + a * x / static_cast<T>(m--);
                           });
  }

  // Returns integral over [a,b].
  T Integrate(T a, T b) const { return Primative(b) - Primative(a); }

 private:
  std::vector<T> coef;
};

}  // namespace Interpolation

#endif  // INTERPOLATION_POLYNOMIAL_GUARD_H
