#ifndef INTERPOLATION_AKIMA_SPLINE_HPP
#define INTERPOLATION_AKIMA_SPLINE_HPP

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief Piecewise-cubic interpolation using locally weighted secant slopes.
 *
 * The object stores iterators rather than copying the samples. The referenced
 * containers must remain alive and must not be reallocated while the
 * interpolator is in use.
 *
 * @tparam xIter Random-access iterator over real abscissae.
 * @tparam yIter Random-access iterator over real or complex ordinates.
 * @pre The abscissa range contains at least three strictly increasing values.
 * Queries outside the sample interval are extrapolated using the first or
 * final cubic segment.
 */
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class Akima {
  public:
    /** @brief Scalar type used for abscissae. */
    using x_value_type = std::iter_value_t<xIter>;
    /** @brief Scalar type used for interpolated values. */
    using y_value_type = std::iter_value_t<yIter>;

    /**
     * @brief Construct an interpolator over non-owning sample ranges.
     * @param xStart Iterator to the first abscissa.
     * @param xFinish Iterator one past the final abscissa.
     * @param yStart Iterator to the ordinate corresponding to `xStart`.
     */
    Akima(xIter xStart, xIter xFinish, yIter yStart);

    /**
     * @brief Evaluate the piecewise-cubic interpolant.
     * @param x Query abscissa.
     * @return Interpolated or extrapolated ordinate.
     */
    y_value_type operator()(x_value_type x) const;

    /**
     * @brief Evaluate the first derivative of the interpolant.
     * @param x Query abscissa.
     * @return First derivative at `x`.
     */
    y_value_type Derivative(x_value_type x) const;

    /**
     * @brief Compatibility spelling for Derivative().
     * @param x Query abscissa.
     * @return First derivative at `x`.
     * @deprecated Use Derivative().
     */
    [[deprecated("Use Derivative()")]] y_value_type deriv(x_value_type x) const {
        return Derivative(x);
    }

  private:
    // Iterators to the function data.
    xIter _xS;
    xIter _xF;
    yIter _yS;

    // m and s values
    std::vector<y_value_type> m;
    std::vector<y_value_type> s;

    std::ptrdiff_t intervalIndex(x_value_type x) const;
};

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::Akima(xIter xS, xIter xF, yIter yS)
    : _xS{xS}, _xF{xF}, _yS{yS} {
    // Dimension of the linear system.
    const auto n = std::distance(_xS, _xF);
    assert(n > 2);
    assert(std::is_sorted(_xS, _xF));

    // fill out m
    m.reserve(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        m.push_back((_yS[i + 1] - _yS[i]) / (_xS[i + 1] - _xS[i]));
    }

    // fill out s
    using weight_type = decltype(std::abs(y_value_type{}));
    constexpr auto oneHalf = static_cast<weight_type>(0.5);
    s.reserve(n);
    s.push_back(m[0]);
    s.push_back((m[0] + m[1]) * oneHalf);
    for (int i = 2; i < n - 2; ++i) {
        const auto a = std::abs(m[i + 1] - m[i]);
        const auto b = std::abs(m[i - 1] - m[i - 2]);
        const auto weightSum = a + b;
        if (weightSum == weight_type{}) {
            s.push_back((m[i - 1] + m[i]) * oneHalf);
        } else {
            s.push_back((a * m[i - 1] + b * m[i]) / weightSum);
        }
    }
    if (n > 3) {
        s.push_back((m[n - 3] + m[n - 2]) * oneHalf);
    }
    s.push_back(m[n - 2]);
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
std::ptrdiff_t
Akima<xIter, yIter>::intervalIndex(x_value_type x) const {
    const auto upper = std::upper_bound(_xS, _xF, x);
    if (upper == _xS) {
        return 0;
    }
    if (upper == _xF) {
        return std::distance(_xS, _xF) - 2;
    }
    return std::distance(_xS, upper) - 1;
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::y_value_type
Akima<xIter, yIter>::operator()(x_value_type x) const {
    const auto i = intervalIndex(x);
    const auto h = _xS[i + 1] - _xS[i];
    const auto offset = x - _xS[i];
    const auto a = _yS[i];
    const auto b = s[i];
    const auto c = (static_cast<y_value_type>(3) * m[i] -
                    static_cast<y_value_type>(2) * s[i] - s[i + 1]) /
                   h;
    const auto d =
        (s[i] + s[i + 1] - static_cast<y_value_type>(2) * m[i]) / (h * h);
    return a + offset * (b + offset * (c + d * offset));
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::y_value_type
Akima<xIter, yIter>::Derivative(x_value_type x) const {
    const auto i = intervalIndex(x);
    const auto h = _xS[i + 1] - _xS[i];
    const auto offset = x - _xS[i];
    const auto b = s[i];
    const auto c = (static_cast<y_value_type>(3) * m[i] -
                    static_cast<y_value_type>(2) * s[i] - s[i + 1]) /
                   h;
    const auto d =
        (s[i] + s[i + 1] - static_cast<y_value_type>(2) * m[i]) / (h * h);
    return b + offset * (static_cast<y_value_type>(2) * c +
                         static_cast<y_value_type>(3) * d * offset);
}

}   // namespace Interpolation

#endif   // INTERPOLATION_AKIMA_SPLINE_HPP
