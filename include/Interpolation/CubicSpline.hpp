#ifndef INTERPOLATION_CUBIC_SPLINE_HPP
#define INTERPOLATION_CUBIC_SPLINE_HPP

#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/CubicSplineSystem.hpp>
#include <Interpolation/Samples.hpp>
#include <Interpolation/Side.hpp>

namespace Interpolation {

/**
 * @brief Piecewise-cubic spline interpolation on ordered sample points.
 *
 * Construction takes ranges, borrowing lvalues and owning rvalues. Queries
 * outside the sample interval continue the first or final cubic piece.
 *
 * The linear system is assembled and factorised by CubicSplineSystem, which
 * this class holds and exposes through System(). A caller fitting many
 * ordinate sets on one grid should build that system once and solve into its
 * own buffers rather than constructing a spline per dataset — the matrix
 * depends on the nodes alone, so there is nothing to recompute.
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
    /** @brief The factorised system this spline was solved with. */
    using System = CubicSplineSystem<XView>;

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
        : _system{Checked(std::move(x), y, left, right), left, right},
          _y{std::move(y)} {
        _ypp.assign(_system.Size(), Scalar{});
        _system.Solve(_y, leftDerivative, rightDerivative,
                      std::span<Scalar>{_ypp});
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
    std::size_t Size() const { return _system.Size(); }

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
        const auto i = Detail::LocateSegment(_system.Nodes(), x);
        return Detail::SplinePiece<N, Real, Scalar>(
            _system.Spacing(i), x - _system.Node(i), _y[i], _y[i + 1], _ypp[i],
            _ypp[i + 1]);
    }

    /** @brief Evaluate the spline; the same as `Evaluate<0>`. */
    Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /**
     * @brief Evaluate the spline, or its `N`th derivative, at every node.
     *
     * Writes one value per node into `out`, in node order, without a segment
     * search: the segment adjoining each node is known. This is what a
     * differentiation operator on a fixed grid wants, and it is where the
     * general query path is doing work it need not.
     *
     * The third derivative jumps across an interior knot; `side` chooses which
     * limit is reported, and defaults to the right-hand one, so that
     * `EvaluateAtNodes<N>(out)` agrees with `Evaluate<N>(Node(k))` for every
     * `k`. Lower orders agree from both sides to rounding.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param out Output, one per node. Overwritten.
     * @param side Which adjoining segment to answer from at a node.
     * @throws std::invalid_argument if `out` has the wrong length.
     */
    template <std::size_t N = 0>
    void EvaluateAtNodes(std::span<Scalar> out, Side side = Side::Right) const {
        _system.template EvaluateAtNodes<N>(_y, _ypp, out, side);
    }

    /**
     * @brief The nodal values or `N`th derivatives, in a fresh vector.
     *
     * The convenient form of EvaluateAtNodes. It allocates, so a caller
     * sweeping many lines should keep one buffer and call EvaluateAtNodes.
     */
    template <std::size_t N = 0>
    std::vector<Scalar> NodeValues(Side side = Side::Right) const {
        std::vector<Scalar> values(Size());
        EvaluateAtNodes<N>(std::span<Scalar>{values}, side);
        return values;
    }

    /** @brief Abscissa of node `i`. */
    Real Node(std::size_t i) const { return _system.Node(i); }

    /** @brief Index of the segment used to evaluate `x`. */
    std::size_t Segment(Real x) const {
        return Detail::LocateSegment(_system.Nodes(), x);
    }

    /**
     * @brief The factorised system behind this spline.
     *
     * Exposed so that a caller who built one spline and then finds it has more
     * data on the same grid can reuse the factorisation rather than pay for it
     * again.
     */
    const System &SplineSystem() const { return _system; }

    /** @brief The nodal second derivatives. */
    std::span<const Scalar> Curvatures() const {
        return std::span<const Scalar>{_ypp};
    }

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
        const auto h = _system.Spacing(i);
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
    System _system;
    YView _y;
    std::vector<Scalar> _ypp; // Second derivatives at the nodes.

    // Validate the sample pair before the system sees the abscissae on their
    // own, so that a length mismatch and a bad grid are both reported against
    // this class rather than against the system it delegates to.
    static XView Checked(XView x, const YView &y, BoundaryCondition left,
                         BoundaryCondition right) {
        const auto notAKnot = left == BoundaryCondition::NotAKnot;
        if (notAKnot != (right == BoundaryCondition::NotAKnot)) {
            throw std::invalid_argument(
                "CubicSpline: NotAKnot constrains the whole system rather "
                "than one endpoint, so it must be used at both ends or "
                "neither");
        }
        Detail::ValidateSamples(x, y, notAKnot ? 4 : 2, "CubicSpline");
        return x;
    }
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y,
          typename... Rest>
CubicSpline(X &&, Y &&,
            Rest...) -> CubicSpline<std::views::all_t<X>, std::views::all_t<Y>>;

} // namespace Interpolation

#endif // INTERPOLATION_CUBIC_SPLINE_HPP
