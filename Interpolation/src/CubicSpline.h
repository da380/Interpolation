#ifndef INTERPOLATION_CUBIC_SPLINE_GUARD_H
#define INTERPOLATION_CUBIC_SPLINE_GUARD_H

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interpolation {

// Enum class for boundary condition types.
enum class CubicSplineBC { Free, Clamped };

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class CubicSpline {
  public:
    // Define some class member types
    using x_value_type = std::iter_value_t<xIter>;
    using y_value_type = std::iter_value_t<yIter>;

    // Default constructor.
    CubicSpline() = default;

    // General constructor.
    CubicSpline(xIter, xIter, yIter, CubicSplineBC, y_value_type, CubicSplineBC,
                y_value_type);

    // Constructor for natural spline.
    CubicSpline(xIter, xIter, yIter);

    // Constructor when boundary conditions are the same.
    CubicSpline(xIter, xIter, yIter, CubicSplineBC, y_value_type, y_value_type);

    // Evaluate interpolating function.
    y_value_type operator()(x_value_type) const;
    y_value_type Derivative(x_value_type) const;

  private:
    using Vector = Eigen::Matrix<y_value_type, Eigen::Dynamic, 1>;
    using Matrix = Eigen::SparseMatrix<y_value_type>;

    xIter _xS;   // Iterator to the start of the x-values.
    xIter _xF;   // Iterator to the end of hte x-values.
    yIter _yS;   // Iterator to the start of the y-values.

    Vector _ypp;   // Cubic spline coefficients
};

// Definition of the main constructor.
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS,
                                       CubicSplineBC left, y_value_type ypl,
                                       CubicSplineBC right, y_value_type ypr)
    : _xS{xS}, _xF{xF}, _yS{yS} {
    // Dimension of the linear system.
    const auto n = std::distance(_xS, _xF);
    assert(std::is_sorted(_xS, _xF));

    // Set up the sparse matrix.
    Matrix A(n, n);
    A.reserve(Eigen::VectorXi::Constant(n, 3));

    // set some constants
    constexpr auto oneThird =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(3);
    constexpr auto oneSixth =
        static_cast<x_value_type>(1) / static_cast<x_value_type>(6);

    // Store only the lower triangle consumed by SimplicialLDLT. A Free
    // endpoint fixes its second derivative to zero, so the adjacent term can
    // be eliminated from the neighbouring equation. Omitting that edge makes
    // the full-size system symmetric without changing its solution.
    for (int i = 0; i < n - 1; ++i) {
        const bool adjacentToFreeLeft = i == 0 && left == CubicSplineBC::Free;
        const bool adjacentToFreeRight =
            i == n - 2 && right == CubicSplineBC::Free;
        if (!adjacentToFreeLeft && !adjacentToFreeRight) {
            A.insert(i + 1, i) = oneSixth * (_xS[i + 1] - _xS[i]);
        }
    }

    // Add in the diagonal.
    if (left == CubicSplineBC::Free) {
        A.insert(0, 0) = 1;
    } else {
        A.insert(0, 0) = oneThird * (_xS[1] - _xS[0]);
    }
    for (int i = 1; i < n - 1; ++i) {
        A.insert(i, i) = oneThird * (_xS[i + 1] - _xS[i - 1]);
    }
    if (right == CubicSplineBC::Free) {
        A.insert(n - 1, n - 1) = 1;
    } else {
        A.insert(n - 1, n - 1) = oneThird * (_xS[n - 1] - _xS[n - 2]);
    }

    // Finalise the matrix construction.
    A.makeCompressed();

    // Set the right hand side.
    Vector rhs(n);
    if (left == CubicSplineBC::Free) {
        rhs(0) = 0;
    } else {
        rhs(0) = (_yS[1] - _yS[0]) / (_xS[1] - _xS[0]) - ypl;
    }
    for (int i = 1; i < n - 1; i++) {
        rhs(i) = (_yS[i + 1] - _yS[i]) / (_xS[i + 1] - _xS[i]) -
                 (_yS[i] - _yS[i - 1]) / (_xS[i] - _xS[i - 1]);
    }
    if (right == CubicSplineBC::Free) {
        rhs(n - 1) = 0;
    } else {
        rhs(n - 1) =
            ypr - (_yS[n - 1] - _yS[n - 2]) / (_xS[n - 1] - _xS[n - 2]);
    }

    // Solve the linear system.
    // A now stores the lower triangle of a symmetric positive-definite system.
    Eigen::SimplicialLDLT<Matrix, Eigen::Lower> solver;
    solver.compute(A);
    _ypp = solver.solve(rhs);
    assert(solver.info() == Eigen::Success);
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
           ((a * a * a - a) * _ypp(i1) + (b * b * b - b) * _ypp(i2)) * h * h *
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
               ((-3 * a * a + 1) * _ypp(i1) + (3 * b * b - 1) * _ypp(i2));
};

}   // namespace Interpolation

#endif   //  INTERPOLATION_CUBIC_SPLINE_GUARD_H
