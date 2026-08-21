#ifndef INTERPOLATION_CUBIC_SPLINE_HPP
#define INTERPOLATION_CUBIC_SPLINE_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <span>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Tridiagonal.hpp>

namespace Interpolation {

/** @brief Boundary-condition types supported by CubicSpline. */
enum class CubicSplineBC {
    /** Natural condition: the endpoint second derivative is zero. */
    Free,
    /** The endpoint first derivative is supplied by the caller. */
    Clamped
};

/**
 * @brief Piecewise-cubic spline interpolation on ordered sample points.
 *
 * The object stores iterators rather than copying the samples. The referenced
 * containers must remain alive and must not be reallocated while the spline is
 * in use. Queries outside the sample interval use the first or final cubic
 * segment for extrapolation.
 *
 * @tparam xIter Random-access iterator over real abscissae.
 * @tparam yIter Random-access iterator over real or complex ordinates.
 * @pre The abscissa range contains at least two strictly increasing values.
 */
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class CubicSpline {
  public:
    /** @brief Scalar type used for abscissae. */
    using x_value_type = std::iter_value_t<xIter>;
    /** @brief Scalar type used for interpolated values. */
    using y_value_type = std::iter_value_t<yIter>;

    /**
     * @brief Construct an unbound spline for later assignment.
     * @warning Evaluation is invalid until a sample-backed spline is assigned.
     */
    CubicSpline() = default;

    /**
     * @brief Construct a spline with independently selected endpoint
     * conditions.
     * @param xStart Iterator to the first abscissa.
     * @param xFinish Iterator one past the final abscissa.
     * @param yStart Iterator to the ordinate corresponding to `xStart`.
     * @param left Left endpoint boundary-condition type.
     * @param leftDerivative Required first derivative when `left` is Clamped;
     *        ignored when it is Free.
     * @param right Right endpoint boundary-condition type.
     * @param rightDerivative Required first derivative when `right` is Clamped;
     *        ignored when it is Free.
     */
    CubicSpline(xIter xStart, xIter xFinish, yIter yStart, CubicSplineBC left,
                y_value_type leftDerivative, CubicSplineBC right,
                y_value_type rightDerivative);

    /**
     * @brief Construct a natural spline with Free conditions at both endpoints.
     * @param xStart Iterator to the first abscissa.
     * @param xFinish Iterator one past the final abscissa.
     * @param yStart Iterator to the ordinate corresponding to `xStart`.
     */
    CubicSpline(xIter xStart, xIter xFinish, yIter yStart);

    /**
     * @brief Construct a spline using one condition type at both endpoints.
     * @param xStart Iterator to the first abscissa.
     * @param xFinish Iterator one past the final abscissa.
     * @param yStart Iterator to the ordinate corresponding to `xStart`.
     * @param both Boundary-condition type for both endpoints.
     * @param leftDerivative Left derivative when `both` is Clamped; otherwise
     *        ignored.
     * @param rightDerivative Right derivative when `both` is Clamped; otherwise
     *        ignored.
     */
    CubicSpline(xIter xStart, xIter xFinish, yIter yStart, CubicSplineBC both,
                y_value_type leftDerivative, y_value_type rightDerivative);

    /**
     * @brief Evaluate the spline interpolant or its endpoint-segment
     * extrapolation.
     * @param x Query abscissa.
     * @return Spline value at `x`.
     */
    y_value_type operator()(x_value_type x) const;

    /**
     * @brief Evaluate the first derivative of the spline.
     * @param x Query abscissa.
     * @return First derivative at `x`.
     */
    y_value_type Derivative(x_value_type x) const;

  private:
    xIter _xS; // Iterator to the start of the x-values.
    xIter _xF; // Iterator to the end of the x-values.
    yIter _yS; // Iterator to the start of the y-values.

    std::vector<y_value_type> _ypp; // Second derivatives at the nodes.
};

// Definition of the main constructor.
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS,
                                       CubicSplineBC left, y_value_type ypl,
                                       CubicSplineBC right, y_value_type ypr)
    : _xS{xS}, _xF{xF}, _yS{yS} {
    const auto n = static_cast<std::size_t>(std::distance(_xS, _xF));
    assert(n >= 2);
    assert(std::is_sorted(_xS, _xF));

    constexpr auto oneThird =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(3);
    constexpr auto oneSixth =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(6);

    // The system for the nodal second derivatives is tridiagonal, so it is
    // assembled as three diagonals and solved directly. Its coefficients are
    // real even when the ordinates are complex, so only the right-hand side
    // carries the ordinate type.
    //
    // The right-hand side is built in _ypp, which the solve then overwrites
    // with the solution.
    std::vector<x_value_type> sub(n, 0), diag(n, 0), super(n, 0);
    _ypp.assign(n, y_value_type{});

    // Interior rows express continuity of the first derivative at each node.
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const auto hPrev = _xS[i] - _xS[i - 1];
        const auto hNext = _xS[i + 1] - _xS[i];
        sub[i] = oneSixth * hPrev;
        diag[i] = oneThird * (hPrev + hNext);
        super[i] = oneSixth * hNext;
        _ypp[i] = (_yS[i + 1] - _yS[i]) / hNext - (_yS[i] - _yS[i - 1]) / hPrev;
    }

    // Left endpoint. A Free condition states directly that the second
    // derivative there is zero.
    if (left == CubicSplineBC::Free) {
        diag[0] = 1;
    } else {
        const auto h = _xS[1] - _xS[0];
        diag[0] = oneThird * h;
        super[0] = oneSixth * h;
        _ypp[0] = (_yS[1] - _yS[0]) / h - ypl;
    }

    // Right endpoint.
    if (right == CubicSplineBC::Free) {
        diag[n - 1] = 1;
    } else {
        const auto h = _xS[n - 1] - _xS[n - 2];
        sub[n - 1] = oneSixth * h;
        diag[n - 1] = oneThird * h;
        _ypp[n - 1] = ypr - (_yS[n - 1] - _yS[n - 2]) / h;
    }

    Detail::SolveTridiagonal<x_value_type, y_value_type>(
        std::span<const x_value_type>{sub}, std::span<x_value_type>{diag},
        std::span<const x_value_type>{super}, std::span<y_value_type>{_ypp});
}

// Definition of the constructor for natural splines.
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS)
    : CubicSpline(xS, xF, yS, CubicSplineBC::Free, 0, CubicSplineBC::Free, 0) {}

// Definition of the constructor when boundary conditions are the same.
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS,
                                       CubicSplineBC both, y_value_type ypl,
                                       y_value_type ypr)
    : CubicSpline(xS, xF, yS, both, ypl, both, ypr) {}

// Evaluation of the interpolating function.
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::y_value_type
CubicSpline<xIter, yIter>::operator()(x_value_type x) const {
    // Find the first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    constexpr auto oneSixth =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(6);
    auto i2 = std::distance(_xS, iter);
    auto i1 = i2 - 1;
    auto x1 = _xS[i1];
    auto x2 = _xS[i2];
    auto h = x2 - x1;
    auto a = (x2 - x) / h;
    auto b = (x - x1) / h;
    return a * _yS[i1] + b * _yS[i2] +
           ((a * a * a - a) * _ypp[i1] + (b * b * b - b) * _ypp[i2]) * h * h *
               oneSixth;
};

// Evaluation of the derivative
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::y_value_type
CubicSpline<xIter, yIter>::Derivative(x_value_type x) const {
    // Find the first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    constexpr auto oneSixth =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(6);
    auto i2 = std::distance(_xS, iter);
    auto i1 = i2 - 1;
    auto x1 = _xS[i1];
    auto x2 = _xS[i2];
    auto h = x2 - x1;
    auto a = (x2 - x) / h;
    auto b = (x - x1) / h;
    return (_yS[i2] - _yS[i1]) / h +
           oneSixth * h *
               ((-3 * a * a + 1) * _ypp[i1] + (3 * b * b - 1) * _ypp[i2]);
};

} // namespace Interpolation

#endif //  INTERPOLATION_CUBIC_SPLINE_HPP
