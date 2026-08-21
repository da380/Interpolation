#ifndef INTERPOLATION_LINEAR_HPP
#define INTERPOLATION_LINEAR_HPP

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <vector>

#include <Interpolation/Concepts.hpp>

namespace Interpolation {

/**
 * @brief Piecewise-linear interpolation over ordered sample points.
 *
 * The object stores iterators rather than copying the samples. The referenced
 * containers must therefore remain alive and must not be reallocated while the
 * interpolator is in use. Queries outside the sample interval are extrapolated
 * using the first or final line segment.
 *
 * @tparam xIter Random-access iterator over real abscissae.
 * @tparam yIter Random-access iterator over real or complex ordinates.
 *
 * @pre The abscissa range contains at least two strictly increasing values.
 */
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class Linear {
  public:
    /** @brief Scalar type used for abscissae. */
    using x_value_type = std::iter_value_t<xIter>;
    /** @brief Scalar type used for interpolated values. */
    using y_value_type = std::iter_value_t<yIter>;

    /**
     * @brief Construct an interpolator over a non-owning sample range.
     * @param xStart Iterator to the first abscissa.
     * @param xFinish Iterator one past the final abscissa.
     * @param yStart Iterator to the ordinate corresponding to `xStart`.
     */
    Linear(xIter xStart, xIter xFinish, yIter yStart);

    /**
     * @brief Evaluate the piecewise-linear interpolant or extrapolant.
     * @param x Query abscissa.
     * @return Interpolated ordinate.
     * @par Complexity
     * Logarithmic search in the number of samples.
     */
    y_value_type operator()(x_value_type x) const;

    /**
     * @brief Evaluate the slope of the selected line segment.
     *
     * At an interior knot the segment to the right is selected; at the final
     * knot the final segment is selected.
     *
     * @param x Query abscissa.
     * @return Piecewise-constant first derivative.
     * @par Complexity
     * Logarithmic search in the number of samples.
     */
    y_value_type Derivative(x_value_type x) const;

  private:
    xIter _xS;   // Iterator to start of x values.
    xIter _xF;   // Iterator to end of x values.
    yIter _yS;   // Iterator to start of y values.
};

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::Linear(xIter xS, xIter xF, yIter yS)
    : _xS{xS}, _xF{xF}, _yS{yS} {}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::y_value_type
Linear<xIter, yIter>::operator()(const x_value_type x) const {
    // Find first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    auto i2 = std::distance(_xS, iter);
    auto i1 = i2 - 1;
    auto x1 = _xS[i1];
    auto x2 = _xS[i2];
    auto h = x2 - x1;
    auto a = (x2 - x) / h;
    auto b = (x - x1) / h;
    return a * _yS[i1] + b * _yS[i2];
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::y_value_type
Linear<xIter, yIter>::Derivative(const x_value_type x) const {
    // Find first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    auto i2 = std::distance(_xS, iter);
    auto i1 = i2 - 1;
    auto x1 = _xS[i1];
    auto x2 = _xS[i2];
    auto h = x2 - x1;
    return (_yS[i2] - _yS[i1]) / h;
}

}   // namespace Interpolation

#endif   //  INTERPOLATION_LINEAR_HPP
