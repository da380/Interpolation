#ifndef INTERPOLATION_LINEAR_HPP
#define INTERPOLATION_LINEAR_HPP

#include <cstddef>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include <Interpolation/Concepts.hpp>
#include <Interpolation/Samples.hpp>
#include <Interpolation/Side.hpp>

namespace Interpolation {

/**
 * @brief Piecewise-linear interpolation over ordered sample points.
 *
 * Construction takes ranges. An lvalue container is borrowed and an rvalue is
 * taken by value, which follows from storing `std::views::all_t`: an lvalue
 * deduces to a `ref_view` and an rvalue to an `owning_view`. There is no
 * policy tag and no second constructor.
 *
 * @code
 * Linear borrowing{x, y};                        // borrows both containers
 * Linear owning{std::move(x), std::move(y)};     // owns both
 * @endcode
 *
 * Queries outside the sample interval continue the first or final segment.
 *
 * @tparam XView View over real abscissae.
 * @tparam YView View over real or complex ordinates.
 */
template <typename XView, typename YView>
    requires InterpolationRanges<XView, YView> && std::ranges::view<XView> &&
             std::ranges::view<YView>
class Linear {
  public:
    /** @brief Abscissa precision. */
    using Real = std::ranges::range_value_t<XView>;
    /** @brief Ordinate type, real or complex. */
    using Scalar = std::ranges::range_value_t<YView>;

    /**
     * @brief Construct from abscissa and ordinate ranges.
     * @param x Strictly increasing abscissae, at least two.
     * @param y Ordinates, the same length as `x`.
     * @throws std::invalid_argument if the ranges differ in length, are too
     *         short, or the abscissae are not strictly increasing.
     */
    constexpr Linear(XView x, YView y) : _x{std::move(x)}, _y{std::move(y)} {
        Detail::ValidateSamples(_x, _y, 2, "Linear");
    }

    /** @brief Number of interpolation nodes. */
    constexpr std::size_t Size() const {
        return static_cast<std::size_t>(std::ranges::size(_x));
    }

    /**
     * @brief Evaluate the interpolant or its `N`th derivative.
     *
     * The interpolant is piecewise linear, so derivatives of order two and
     * above are identically zero. At an interior knot the segment to the
     * right is used.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param x Query abscissa.
     */
    template <std::size_t N = 0> constexpr Scalar Evaluate(Real x) const {
        const auto i = Detail::LocateSegment(_x, x);
        return Detail::LinearPiece<N, Real, Scalar>(
            _x[i + 1] - _x[i], x - _x[i], _y[i], _y[i + 1]);
    }

    /** @brief Evaluate the interpolant; the same as `Evaluate<0>`. */
    constexpr Scalar operator()(Real x) const { return Evaluate<0>(x); }

    /**
     * @brief Evaluate the interpolant, or its `N`th derivative, at every node.
     *
     * Writes one value per node into `out`, in node order, without a segment
     * search: the segment adjoining each node is known. For a caller who wants
     * the whole nodal sweep — a difference operator on a fixed grid, say —
     * this replaces `Size()` binary searches with none.
     *
     * The first derivative jumps across an interior node, since the pieces are
     * linear; `side` chooses which limit is reported, and defaults to the
     * right-hand one, so that `EvaluateAtNodes<N>(out)` agrees with
     * `Evaluate<N>(Node(k))` for every `k`.
     *
     * @tparam N Derivative order; `0` is the value itself.
     * @param out Output, one per node. Overwritten.
     * @param side Which adjoining segment to answer from at a node.
     * @throws std::invalid_argument if `out` has the wrong length.
     */
    template <std::size_t N = 0>
    void EvaluateAtNodes(std::span<Scalar> out, Side side = Side::Right) const {
        const auto n = Size();
        Detail::ValidateNodeOutput(out.size(), n, "Linear");
        for (std::size_t k = 0; k < n; ++k) {
            const auto i = Detail::NodeSegment(n, k, side);
            out[k] = Detail::LinearPiece<N, Real, Scalar>(
                _x[i + 1] - _x[i], _x[k] - _x[i], _y[i], _y[i + 1]);
        }
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
    constexpr Real Node(std::size_t i) const { return _x[i]; }

    /** @brief Index of the segment used to evaluate `x`. */
    constexpr std::size_t Segment(Real x) const {
        return Detail::LocateSegment(_x, x);
    }

    /**
     * @brief Integral of segment `i` from its left node over a width `t`.
     *
     * Exposed so that Primitive() can accumulate an antiderivative without
     * this class having to carry one. On a segment the interpolant is
     * `y[i] + m t`, so the integral is `y[i] t + m t^2 / 2`.
     */
    constexpr Scalar SegmentIntegral(std::size_t i, Real t) const {
        const auto h = _x[i + 1] - _x[i];
        const auto slope = (_y[i + 1] - _y[i]) / h;
        return _y[i] * t + slope * (t * t) / static_cast<Real>(2);
    }

  private:
    XView _x;
    YView _y;
};

/// Borrow lvalue containers and own rvalue ones.
template <std::ranges::viewable_range X, std::ranges::viewable_range Y>
Linear(X &&, Y &&) -> Linear<std::views::all_t<X>, std::views::all_t<Y>>;

} // namespace Interpolation

#endif // INTERPOLATION_LINEAR_HPP
