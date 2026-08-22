#ifndef INTERPOLATION_CUBIC_SPLINE_HPP
#define INTERPOLATION_CUBIC_SPLINE_HPP

#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Samples.hpp>
#include <Interpolation/Tridiagonal.hpp>

namespace Interpolation {

/** @brief Endpoint conditions for CubicSpline. */
enum class BoundaryCondition {
    /** The endpoint second derivative is zero. */
    Natural,
    /** The endpoint first derivative is supplied by the caller. */
    Clamped,
    /**
     * @brief The third derivative is continuous across the first and last
     * interior knots.
     *
     * Imposes nothing false at the boundary, so a cubic is reproduced exactly
     * and the scheme stays fourth order right up to the ends, where Natural
     * costs an order. It is a condition on the whole system rather than on one
     * endpoint, so it must be used at both ends and needs at least four nodes.
     */
    NotAKnot
};

/**
 * @brief Piecewise-cubic spline interpolation on ordered sample points.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues. Queries
 * outside the sample interval continue the first or final cubic piece.
 *
 * @tparam XView View over real abscissae.
 * @tparam YView View over real or complex ordinates.
 */
template <typename XView, typename YView>
    requires InterpolationRanges<XView, YView> && std::ranges::view<XView> &&
             std::ranges::view<YView>
class CubicSpline {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;
    /** @brief Ordinate type, real or complex. */
    using Scalar = std::ranges::range_value_t<YView>;

    /**
     * @brief Construct a spline with independently chosen endpoint conditions.
     * @param x Strictly increasing abscissae, at least two.
     * @param y Ordinates, the same length as `x`.
     * @param left Condition at the first node.
     * @param leftDerivative First derivative there when `left` is Clamped.
     * @param right Condition at the final node.
     * @param rightDerivative First derivative there when `right` is Clamped.
     * @throws std::invalid_argument if the ranges differ in length, are too
     *         short, or the abscissae are not strictly increasing.
     */
    CubicSpline(XView x, YView y, BoundaryCondition left, Scalar leftDerivative,
                BoundaryCondition right, Scalar rightDerivative)
        : _x{std::move(x)}, _y{std::move(y)} {
        const auto notAKnot = left == BoundaryCondition::NotAKnot;
        if (notAKnot != (right == BoundaryCondition::NotAKnot)) {
            throw std::invalid_argument(
                "CubicSpline: NotAKnot constrains the whole system rather "
                "than one endpoint, so it must be used at both ends or "
                "neither");
        }
        Detail::ValidateSamples(_x, _y, notAKnot ? 4 : 2, "CubicSpline");
        if (notAKnot) {
            SolveNotAKnot();
        } else {
            Solve(left, leftDerivative, right, rightDerivative);
        }
    }

    /** @brief Construct a natural spline, with Natural at both endpoints. */
    CubicSpline(XView x, YView y)
        : CubicSpline{std::move(x),
                      std::move(y),
                      BoundaryCondition::Natural,
                      Scalar{},
                      BoundaryCondition::Natural,
                      Scalar{}} {}

    /** @brief Construct with the same condition type at both endpoints. */
    CubicSpline(XView x, YView y, BoundaryCondition both, Scalar leftDerivative,
                Scalar rightDerivative)
        : CubicSpline{std::move(x),   std::move(y), both,
                      leftDerivative, both,         rightDerivative} {}

    /** @brief Number of interpolation nodes. */
    std::size_t Size() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /**
     * @brief Evaluate the spline or its `N`th derivative.
     *
     * The pieces are cubic, so derivatives of order four and above are
     * identically zero and the third is piecewise constant.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> Scalar Evaluate(Real x) const {
        const auto i = Detail::LocateSegment(_x, x);
        return Detail::SplinePiece<N, Real, Scalar>(_x[i + 1] - _x[i],
                                                    x - _x[i], _y[i], _y[i + 1],
                                                    _ypp[i], _ypp[i + 1]);
    }

    /** @brief Evaluate the spline; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /** @brief Abscissa of node `i`. */
    Real Node(std::size_t i) const { return _x[i]; }

    /** @brief Index of the segment used to evaluate `x`. */
    std::size_t Segment(Real x) const { return Detail::LocateSegment(_x, x); }

    /**
     * @brief Integral of segment `i` from its left node over a width `t`.
     *
     * With @f$ u = t/h @f$ the segment integrates term by term to
     *
     * @f[
     * h\left[ y_i\left(u - \frac{u^2}{2}\right)
     *   + y_{i+1}\frac{u^2}{2}
     *   + \frac{h^2}{6}\left(
     *       M_i\left(-\frac14 - \frac{(1-u)^4}{4} + \frac{(1-u)^2}{2}\right)
     *     + M_{i+1}\left(\frac{u^4}{4} - \frac{u^2}{2}\right)
     *     \right)\right],
     * @f]
     *
     * which at @f$ u = 1 @f$ reduces to the familiar
     * @f$ h(y_i + y_{i+1})/2 - h^3(M_i + M_{i+1})/24 @f$.
     */
    Scalar SegmentIntegral(std::size_t i, Real t) const {
        const auto h = _x[i + 1] - _x[i];
        const auto u = t / h;
        const auto v = static_cast<Real>(1) - u;
        constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
        constexpr auto quarter = static_cast<Real>(1) / static_cast<Real>(4);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

        const auto left = u - half * u * u;
        const auto right = half * u * u;
        const auto curvatureLeft =
            -quarter - quarter * v * v * v * v + half * v * v;
        const auto curvatureRight = quarter * u * u * u * u - half * u * u;

        return h *
               (_y[i] * left + _y[i + 1] * right +
                oneSixth * h * h *
                    (_ypp[i] * curvatureLeft + _ypp[i + 1] * curvatureRight));
    }

  private:
    XView _x;
    YView _y;
    std::vector<Scalar> _ypp; // Second derivatives at the nodes.

    // Assemble and solve the tridiagonal system for the nodal second
    // derivatives. The coefficients are real even when the ordinates are
    // complex, so only the right-hand side carries the ordinate type. The
    // right-hand side is built in _ypp, which the solve overwrites with the
    // solution.
    // Not-a-knot is solved in slope form and converted, for the reasons set
    // out on Detail::NotAKnotCurvatures.
    void SolveNotAKnot() {
        const auto n = Size();
        std::vector<Real> nodes(n);
        std::vector<Scalar> values(n);
        for (std::size_t i = 0; i < n; ++i) {
            nodes[i] = _x[i];
            values[i] = _y[i];
        }

        _ypp.assign(n, Scalar{});
        std::vector<Real> sub(n), diag(n), super(n);
        std::vector<Scalar> slope(n);

        Detail::NotAKnotCurvatures<Real, Scalar>(
            std::span<const Real>{nodes}, std::span<const Scalar>{values},
            std::span<Scalar>{_ypp}, std::span<Real>{sub},
            std::span<Real>{diag}, std::span<Real>{super},
            std::span<Scalar>{slope});
    }

    void Solve(BoundaryCondition left, Scalar leftDerivative,
               BoundaryCondition right, Scalar rightDerivative) {
        const auto n = Size();
        constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

        std::vector<Real> sub(n, 0), diag(n, 0), super(n, 0);
        _ypp.assign(n, Scalar{});

        // Interior rows express continuity of the first derivative.
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const auto hPrev = _x[i] - _x[i - 1];
            const auto hNext = _x[i + 1] - _x[i];
            sub[i] = oneSixth * hPrev;
            diag[i] = oneThird * (hPrev + hNext);
            super[i] = oneSixth * hNext;
            _ypp[i] = (_y[i + 1] - _y[i]) / hNext - (_y[i] - _y[i - 1]) / hPrev;
        }

        // A Natural condition states directly that the second derivative
        // at that endpoint is zero.
        if (left == BoundaryCondition::Natural) {
            diag[0] = 1;
        } else {
            const auto h = _x[1] - _x[0];
            diag[0] = oneThird * h;
            super[0] = oneSixth * h;
            _ypp[0] = (_y[1] - _y[0]) / h - leftDerivative;
        }

        if (right == BoundaryCondition::Natural) {
            diag[n - 1] = 1;
        } else {
            const auto h = _x[n - 1] - _x[n - 2];
            sub[n - 1] = oneSixth * h;
            diag[n - 1] = oneThird * h;
            _ypp[n - 1] = rightDerivative - (_y[n - 1] - _y[n - 2]) / h;
        }

        Detail::SolveTridiagonal<Real, Scalar>(
            std::span<const Real>{sub}, std::span<Real>{diag},
            std::span<const Real>{super}, std::span<Scalar>{_ypp});
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y,
          typename... Rest>
CubicSpline(X &&, Y &&,
            Rest...) -> CubicSpline<std::views::all_t<X>, std::views::all_t<Y>>;

} // namespace Interpolation

#endif // INTERPOLATION_CUBIC_SPLINE_HPP
